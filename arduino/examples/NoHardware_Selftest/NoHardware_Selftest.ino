/*
 * NoHardware_Selftest — verifies the ArbiterAHRS library without any sensor.
 *
 * The sketch feeds the filter synthetic measurements of a body resting at
 * roll = 10, pitch = 5, yaw = 60 degrees, with a deliberate 30-degree
 * initial heading error and a gyro bias of about 1 deg/s. Expected result
 * after 20 simulated seconds: total attitude error below 1.5 degrees and
 * a bias estimate within 0.15 deg/s of the injected value.
 *
 * No wiring is required. Open the Serial Monitor at 115200 baud.
 */
#include <ArbiterAHRS.h>

static const float DEG = 0.017453293f;
static const float DT  = 0.004f;          /* 250 Hz */

ArbiterAHRS ahrs;

/* True orientation: roll 10, pitch 5, yaw 60 deg (ZYX). */
static void trueVectors(float v_e[3], float out_b[3],
                        float roll, float pitch, float yaw) {
    /* body vector = R^T * earth vector, R built from ZYX Euler angles */
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

static float frand() { return (float)random(-1000, 1001) / 1000.0f; }

void setup() {
    Serial.begin(115200);
    while (!Serial) {}
    randomSeed(12345);
    ahrs.begin();

    const float roll0 = 10*DEG, pitch0 = 5*DEG, yaw0 = 60*DEG;
    float g_e[3] = { 0.0f, 0.0f, 9.81f };
    float m_e[3] = { 0.0f, 0.5f, -0.866f };          /* dip 60 deg */
    float acc_b[3], mag_b[3], mag_bad[3];
    trueVectors(g_e, acc_b, roll0, pitch0, yaw0);
    trueVectors(m_e, mag_b, roll0, pitch0, yaw0);
    /* first magnetometer sample rotated 30 deg about the vertical:
       forces an initial heading error the filter must remove */
    trueVectors(m_e, mag_bad, roll0, pitch0, yaw0 + 30*DEG);

    const float bias[3] = { 1.2f*DEG, -0.8f*DEG, 0.6f*DEG };

    ahrs.update(0, 0, 0, acc_b[0], acc_b[1], acc_b[2],
                mag_bad[0], mag_bad[1], mag_bad[2], DT);

    for (long i = 0; i < (long)(20.0f / DT); i++) {
        float gx = bias[0] + frand()*0.15f*DEG;
        float gy = bias[1] + frand()*0.15f*DEG;
        float gz = bias[2] + frand()*0.15f*DEG;
        float ax = acc_b[0] + frand()*0.04f;
        float ay = acc_b[1] + frand()*0.04f;
        float az = acc_b[2] + frand()*0.04f;
        float mx = mag_b[0] + frand()*0.01f;
        float my = mag_b[1] + frand()*0.01f;
        float mz = mag_b[2] + frand()*0.01f;
        ahrs.update(gx, gy, gz, ax, ay, az, mx, my, mz, DT);
    }

    float r, p, y, bx, by, bz;
    ahrs.getEulerDeg(r, p, y);
    ahrs.getGyroBias(bx, by, bz);
    Serial.println(F("=== ArbiterAHRS self-test ==="));
    Serial.print(F("euler deg (expect ~10, ~5, ~60): "));
    Serial.print(r, 2); Serial.print(' ');
    Serial.print(p, 2); Serial.print(' ');
    Serial.println(y, 2);
    Serial.print(F("gyro bias deg/s (expect ~1.2, ~-0.8, ~0.6): "));
    Serial.print(bx / DEG, 3); Serial.print(' ');
    Serial.print(by / DEG, 3); Serial.print(' ');
    Serial.println(bz / DEG, 3);
    Serial.print(F("accel evidence: ")); Serial.println(ahrs.accelEvidence(), 2);
    Serial.print(F("mag evidence:   ")); Serial.println(ahrs.magEvidence(), 2);

    bool ok = fabsf(r - 10) < 1.5f && fabsf(p - 5) < 1.5f &&
              fabsf(y - 60) < 1.5f && fabsf(bx/DEG - 1.2f) < 0.15f;
    Serial.println(ok ? F("RESULT: PASS") : F("RESULT: FAIL"));
}

void loop() {}
