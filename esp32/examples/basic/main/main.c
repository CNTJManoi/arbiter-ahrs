/*
 * ARBITER on ESP32: self-test without hardware, then (optionally) live
 * fusion with an MPU-9250 over I2C.
 *
 * The self-test feeds synthetic static measurements with a known attitude
 * (roll 10, pitch 5, yaw 60 deg), a 30-degree initial heading error, and a
 * gyro bias of about 1 deg/s. Expected outcome: attitude error below
 * 1.5 degrees and bias recovered within 0.15 deg/s after 20 simulated
 * seconds. It runs on any ESP32 with no wiring.
 *
 * To run the live part, wire an MPU-9250 (SDA GPIO21, SCL GPIO22 by
 * default), set RUN_LIVE_MPU9250 to 1, and rebuild. The magnetometer must
 * be hard-iron calibrated for meaningful heading; subtract your offsets
 * where marked.
 */
#include <math.h>
#include <stdio.h>
#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "esp_timer.h"
#include "esp_log.h"
#include "arbiter.h"

#define RUN_LIVE_MPU9250 0

static const char *TAG = "arbiter";
#define DEG 0.017453293f

/* ------------------------------------------------------------ self-test */
static uint32_t rng_state = 12345;
static float frand(void) {              /* uniform in [-1, 1] */
    rng_state = rng_state * 1664525u + 1013904223u;
    return ((int32_t)rng_state >> 8) / 8388608.0f;
}

static void body_from_earth(const float v_e[3], float out_b[3],
                            float roll, float pitch, float yaw) {
    float cr = cosf(roll), sr = sinf(roll);
    float cp = cosf(pitch), sp = sinf(pitch);
    float cy = cosf(yaw), sy = sinf(yaw);
    float R[3][3] = {
        { cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr },
        { sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr },
        { -sp,   cp*sr,            cp*cr            }};
    for (int i = 0; i < 3; i++)
        out_b[i] = R[0][i]*v_e[0] + R[1][i]*v_e[1] + R[2][i]*v_e[2];
}

static void run_selftest(void) {
    ArbiterParams p;
    ArbiterState s;
    arbiter_default_params(&p);
    arbiter_init(&s, &p);

    const float dt = 0.004f;
    const float roll0 = 10*DEG, pitch0 = 5*DEG, yaw0 = 60*DEG;
    const float g_e[3] = { 0, 0, 9.81f };
    const float m_e[3] = { 0, 0.5f, -0.866f };
    const float bias[3] = { 1.2f*DEG, -0.8f*DEG, 0.6f*DEG };
    float acc_b[3], mag_b[3], mag_bad[3];
    body_from_earth(g_e, acc_b, roll0, pitch0, yaw0);
    body_from_earth(m_e, mag_b, roll0, pitch0, yaw0);
    body_from_earth(m_e, mag_bad, roll0, pitch0, yaw0 + 30*DEG);

    int64_t t0 = esp_timer_get_time();
    arbiter_update(&s, 0, 0, 0, acc_b[0], acc_b[1], acc_b[2],
                   mag_bad[0], mag_bad[1], mag_bad[2], dt);
    long steps = (long)(20.0f / dt);
    for (long i = 0; i < steps; i++) {
        arbiter_update(&s,
            bias[0] + frand()*0.15f*DEG,
            bias[1] + frand()*0.15f*DEG,
            bias[2] + frand()*0.15f*DEG,
            acc_b[0] + frand()*0.04f,
            acc_b[1] + frand()*0.04f,
            acc_b[2] + frand()*0.04f,
            mag_b[0] + frand()*0.01f,
            mag_b[1] + frand()*0.01f,
            mag_b[2] + frand()*0.01f,
            dt);
    }
    int64_t us = esp_timer_get_time() - t0;

    float r, pch, y;
    arbiter_euler(&s, &r, &pch, &y);
    ESP_LOGI(TAG, "self-test euler deg: %.2f %.2f %.2f (expect ~10 ~5 ~60)",
             r/DEG, pch/DEG, y/DEG);
    ESP_LOGI(TAG, "bias deg/s: %.3f %.3f %.3f (expect ~1.2 ~-0.8 ~0.6)",
             s.b[0]/DEG, s.b[1]/DEG, s.b[2]/DEG);
    ESP_LOGI(TAG, "per-step time on this chip: %.2f us",
             (double)us / (double)(steps + 1));
    int ok = fabsf(r/DEG - 10) < 1.5f && fabsf(pch/DEG - 5) < 1.5f &&
             fabsf(y/DEG - 60) < 1.5f;
    ESP_LOGI(TAG, "RESULT: %s", ok ? "PASS" : "FAIL");
}

/* ------------------------------------------------- live MPU-9250 (option) */
#if RUN_LIVE_MPU9250
#include "driver/i2c.h"

#define I2C_PORT I2C_NUM_0
#define PIN_SDA 21
#define PIN_SCL 22
#define MPU 0x68
#define MAG 0x0C
static const float GYRO_SCALE = 500.0f/32768.0f*DEG;
static const float ACC_SCALE  = 4.0f*9.80665f/32768.0f;
static float mag_adj[3];

static void wr(uint8_t dev, uint8_t reg, uint8_t val) {
    uint8_t buf[2] = { reg, val };
    i2c_master_write_to_device(I2C_PORT, dev, buf, 2, pdMS_TO_TICKS(20));
}
static void rd(uint8_t dev, uint8_t reg, uint8_t n, uint8_t *buf) {
    i2c_master_write_read_device(I2C_PORT, dev, &reg, 1, buf, n,
                                 pdMS_TO_TICKS(20));
}

static void mpu9250_init(void) {
    i2c_config_t c = {
        .mode = I2C_MODE_MASTER,
        .sda_io_num = PIN_SDA, .scl_io_num = PIN_SCL,
        .sda_pullup_en = GPIO_PULLUP_ENABLE,
        .scl_pullup_en = GPIO_PULLUP_ENABLE,
        .master.clk_speed = 400000,
    };
    i2c_param_config(I2C_PORT, &c);
    i2c_driver_install(I2C_PORT, I2C_MODE_MASTER, 0, 0, 0);
    wr(MPU, 0x6B, 0x01); vTaskDelay(pdMS_TO_TICKS(50));
    wr(MPU, 0x1A, 0x03); wr(MPU, 0x1B, 0x08);
    wr(MPU, 0x1C, 0x08); wr(MPU, 0x1D, 0x03);
    wr(MPU, 0x37, 0x02);
    uint8_t asa[3];
    wr(MAG, 0x0A, 0x00); vTaskDelay(pdMS_TO_TICKS(10));
    wr(MAG, 0x0A, 0x0F); vTaskDelay(pdMS_TO_TICKS(10));
    rd(MAG, 0x10, 3, asa);
    for (int i = 0; i < 3; i++)
        mag_adj[i] = ((float)asa[i] - 128.0f)/256.0f + 1.0f;
    wr(MAG, 0x0A, 0x00); vTaskDelay(pdMS_TO_TICKS(10));
    wr(MAG, 0x0A, 0x16);
}

static void live_task(void *arg) {
    ArbiterParams p; ArbiterState s;
    arbiter_default_params(&p);
    arbiter_init(&s, &p);
    mpu9250_init();
    int64_t last = esp_timer_get_time();
    int n = 0;
    for (;;) {
        uint8_t b[14], st, mb[7];
        rd(MPU, 0x3B, 14, b);
        int16_t axr = (b[0]<<8)|b[1], ayr = (b[2]<<8)|b[3], azr = (b[4]<<8)|b[5];
        int16_t gxr = (b[8]<<8)|b[9], gyr = (b[10]<<8)|b[11], gzr = (b[12]<<8)|b[13];
        float mx = 0, my = 0, mz = 0;
        rd(MAG, 0x02, 1, &st);
        if (st & 0x01) {
            rd(MAG, 0x03, 7, mb);
            if (!(mb[6] & 0x08)) {
                int16_t x = (mb[1]<<8)|mb[0], y = (mb[3]<<8)|mb[2],
                        z = (mb[5]<<8)|mb[4];
                mx = (float)y * mag_adj[1];   /* TODO: hard-iron offsets */
                my = (float)x * mag_adj[0];
                mz = -(float)z * mag_adj[2];
            }
        }
        int64_t now = esp_timer_get_time();
        float dt = (now - last) * 1e-6f;
        last = now;
        arbiter_update(&s, gxr*GYRO_SCALE, gyr*GYRO_SCALE, gzr*GYRO_SCALE,
                       axr*ACC_SCALE, ayr*ACC_SCALE, azr*ACC_SCALE,
                       mx, my, mz, dt);
        if (++n % 50 == 0) {
            float r, pch, yw;
            arbiter_euler(&s, &r, &pch, &yw);
            ESP_LOGI(TAG, "rpy: %6.1f %6.1f %6.1f  [%c%c]",
                     r/DEG, pch/DEG, yw/DEG,
                     s.gate_a ? 'A' : '-', s.gate_m ? 'M' : '-');
        }
        vTaskDelay(1);
    }
}
#endif /* RUN_LIVE_MPU9250 */

void app_main(void) {
    run_selftest();
#if RUN_LIVE_MPU9250
    xTaskCreate(live_task, "arbiter_live", 4096, NULL, 5, NULL);
#endif
}
