/*
 * selftest_main.c — ARBITER self-test on STM32 without any sensor.
 *
 * Add arbiter.c/arbiter.h to your CubeMX project (Core/Src, Core/Inc or a
 * separate library folder), enable a UART for printf (or watch the
 * variables in the debugger), and call arbiter_selftest() from main()
 * after HAL_Init()/clock setup.
 *
 * The test feeds synthetic static measurements with a known attitude
 * (roll 10, pitch 5, yaw 60 deg), a 30-degree initial heading error, and
 * a gyro bias of about 1 deg/s. Expected outcome after 20 simulated
 * seconds: attitude error below 1.5 deg, bias recovered within
 * 0.15 deg/s. The DWT cycle counter measures the real per-step cost of
 * the filter on your MCU.
 */
#include <math.h>
#include <stdio.h>
#include "arbiter.h"
/* #include "main.h"  -- your CubeMX header, for HAL types if needed */

#define DEG 0.017453293f

static uint32_t rng_state = 12345u;
static float frand(void) {
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

/* Cortex-M3/M4/M7: DWT cycle counter for per-step timing. */
static void dwt_enable(void) {
#ifdef DWT
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->CYCCNT = 0;
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
#endif
}

void arbiter_selftest(void) {
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
    body_from_earth(m_e, mag_b, roll0, pitch0, yaw0 + 0*DEG);
    body_from_earth(m_e, mag_bad, roll0, pitch0, yaw0 + 30*DEG);

    dwt_enable();
    uint32_t cycles = 0;
    arbiter_update(&s, 0, 0, 0, acc_b[0], acc_b[1], acc_b[2],
                   mag_bad[0], mag_bad[1], mag_bad[2], dt);
    long steps = (long)(20.0f / dt);
    for (long i = 0; i < steps; i++) {
        float gx = bias[0] + frand()*0.15f*DEG;
        float gy = bias[1] + frand()*0.15f*DEG;
        float gz = bias[2] + frand()*0.15f*DEG;
        float ax = acc_b[0] + frand()*0.04f;
        float ay = acc_b[1] + frand()*0.04f;
        float az = acc_b[2] + frand()*0.04f;
        float mx = mag_b[0] + frand()*0.01f;
        float my = mag_b[1] + frand()*0.01f;
        float mz = mag_b[2] + frand()*0.01f;
#ifdef DWT
        uint32_t c0 = DWT->CYCCNT;
#endif
        arbiter_update(&s, gx, gy, gz, ax, ay, az, mx, my, mz, dt);
#ifdef DWT
        cycles += DWT->CYCCNT - c0;
#endif
    }

    float r, pch, y;
    arbiter_euler(&s, &r, &pch, &y);
    /* Route printf to your UART (retarget _write) or inspect in debugger. */
    printf("euler deg: %.2f %.2f %.2f (expect ~10 ~5 ~60)\r\n",
           (double)(r/DEG), (double)(pch/DEG), (double)(y/DEG));
    printf("bias deg/s: %.3f %.3f %.3f\r\n",
           (double)(s.b[0]/DEG), (double)(s.b[1]/DEG), (double)(s.b[2]/DEG));
    if (cycles)
        printf("avg cycles/step: %lu\r\n", (unsigned long)(cycles / steps));
    printf("RESULT: %s\r\n",
           (fabsf(r/DEG - 10) < 1.5f && fabsf(pch/DEG - 5) < 1.5f &&
            fabsf(y/DEG - 60) < 1.5f) ? "PASS" : "FAIL");
}
