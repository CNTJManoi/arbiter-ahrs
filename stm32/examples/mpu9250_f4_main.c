/*
 * mpu9250_f4_main.c — ARBITER with an MPU-9250 on STM32F4 (HAL, I2C).
 *
 * Prerequisites (CubeMX): I2C1 at 400 kHz on the pins of your board, a
 * UART for output, and a 1 kHz SysTick (default). Add arbiter.c/arbiter.h
 * to the project. Replace the main() skeleton markers with the CubeMX
 * generated init calls of your project.
 *
 * Sensor configuration: gyro +/-500 dps, accel +/-4 g, AK8963 16-bit at
 * 100 Hz through the bypass mux. The magnetometer axes are remapped into
 * the gyro/accel body frame (x<->y swap, z negated). Hard-iron offsets
 * must be subtracted where marked; without them the heading converges to
 * a distorted north.
 *
 * Loop rate: the code runs the filter on every gyro/accel read (about
 * 500 Hz with a 41 Hz DLPF); dt is measured with the DWT cycle counter,
 * so sampling jitter does not degrade the integration.
 */
#include <math.h>
#include <stdio.h>
#include "main.h"          /* CubeMX: HAL, hi2c1, huart2 */
#include "arbiter.h"

#define DEG 0.017453293f
#define MPU (0x68 << 1)    /* HAL uses 8-bit addresses */
#define MAG (0x0C << 1)

extern I2C_HandleTypeDef hi2c1;
extern UART_HandleTypeDef huart2;

static const float GYRO_SCALE = 500.0f / 32768.0f * DEG;
static const float ACC_SCALE  = 4.0f * 9.80665f / 32768.0f;
static float mag_adj[3];

static void wr(uint16_t dev, uint8_t reg, uint8_t val) {
    HAL_I2C_Mem_Write(&hi2c1, dev, reg, 1, &val, 1, 20);
}
static void rd(uint16_t dev, uint8_t reg, uint8_t n, uint8_t *buf) {
    HAL_I2C_Mem_Read(&hi2c1, dev, reg, 1, buf, n, 20);
}

static void mpu9250_init(void) {
    wr(MPU, 0x6B, 0x01); HAL_Delay(50);   /* wake, PLL clock */
    wr(MPU, 0x1A, 0x03);                  /* gyro DLPF 41 Hz */
    wr(MPU, 0x1B, 0x08);                  /* +/-500 dps */
    wr(MPU, 0x1C, 0x08);                  /* +/-4 g */
    wr(MPU, 0x1D, 0x03);                  /* accel DLPF 41 Hz */
    wr(MPU, 0x37, 0x02);                  /* bypass mux to AK8963 */
    uint8_t asa[3];
    wr(MAG, 0x0A, 0x00); HAL_Delay(10);
    wr(MAG, 0x0A, 0x0F); HAL_Delay(10);
    rd(MAG, 0x10, 3, asa);
    for (int i = 0; i < 3; i++)
        mag_adj[i] = ((float)asa[i] - 128.0f) / 256.0f + 1.0f;
    wr(MAG, 0x0A, 0x00); HAL_Delay(10);
    wr(MAG, 0x0A, 0x16);                  /* 16-bit continuous 100 Hz */
}

static int read_mag(float m[3]) {
    uint8_t st, b[7];
    rd(MAG, 0x02, 1, &st);
    if (!(st & 0x01)) return 0;
    rd(MAG, 0x03, 7, b);
    if (b[6] & 0x08) return 0;            /* magnetic overflow */
    int16_t x = (int16_t)(b[1] << 8 | b[0]);
    int16_t y = (int16_t)(b[3] << 8 | b[2]);
    int16_t z = (int16_t)(b[5] << 8 | b[4]);
    m[0] = (float)y * mag_adj[1];         /* TODO: subtract hard-iron */
    m[1] = (float)x * mag_adj[0];
    m[2] = -(float)z * mag_adj[2];
    return 1;
}

int main(void) {
    HAL_Init();
    /* SystemClock_Config();  -- CubeMX generated */
    /* MX_GPIO_Init(); MX_I2C1_Init(); MX_USART2_UART_Init(); */

    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->CYCCNT = 0;
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;

    ArbiterParams params;
    ArbiterState  ahrs;
    arbiter_default_params(&params);
    arbiter_init(&ahrs, &params);
    mpu9250_init();

    uint32_t last_cyc = DWT->CYCCNT, last_print = HAL_GetTick();
    const float cyc2s = 1.0f / (float)SystemCoreClock;

    for (;;) {
        uint8_t b[14];
        rd(MPU, 0x3B, 14, b);
        int16_t axr = (int16_t)(b[0] << 8 | b[1]);
        int16_t ayr = (int16_t)(b[2] << 8 | b[3]);
        int16_t azr = (int16_t)(b[4] << 8 | b[5]);
        int16_t gxr = (int16_t)(b[8] << 8 | b[9]);
        int16_t gyr = (int16_t)(b[10] << 8 | b[11]);
        int16_t gzr = (int16_t)(b[12] << 8 | b[13]);

        float m[3] = { 0.0f, 0.0f, 0.0f };
        read_mag(m);                       /* zeros mean “no sample” */

        uint32_t now_cyc = DWT->CYCCNT;
        float dt = (float)(now_cyc - last_cyc) * cyc2s;
        last_cyc = now_cyc;
        if (dt <= 0.0f || dt > 0.5f) continue;

        arbiter_update(&ahrs,
                       gxr * GYRO_SCALE, gyr * GYRO_SCALE, gzr * GYRO_SCALE,
                       axr * ACC_SCALE,  ayr * ACC_SCALE,  azr * ACC_SCALE,
                       m[0], m[1], m[2], dt);

        if (HAL_GetTick() - last_print >= 50) {
            last_print = HAL_GetTick();
            float roll, pitch, yaw;
            arbiter_euler(&ahrs, &roll, &pitch, &yaw);
            char line[96];
            int n = snprintf(line, sizeof line,
                             "%6.1f %6.1f %6.1f  [%c%c]\r\n",
                             (double)(roll / DEG), (double)(pitch / DEG),
                             (double)(yaw / DEG),
                             ahrs.gate_a ? 'A' : '-',
                             ahrs.gate_m ? 'M' : '-');
            HAL_UART_Transmit(&huart2, (uint8_t *)line, n, 20);
        }
    }
}
