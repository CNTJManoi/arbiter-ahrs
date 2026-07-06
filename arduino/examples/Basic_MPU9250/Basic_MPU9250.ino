/*
 * Basic_MPU9250 — ArbiterAHRS with an MPU-9250 over I2C (Wire).
 *
 * Wiring: SDA/SCL to the board I2C pins, VCC 3.3 V, AD0 to GND
 * (I2C address 0x68). The AK8963 magnetometer is accessed through
 * the MPU-9250 bypass mux at address 0x0C.
 *
 * Ranges used: gyro +/-500 dps, accel +/-4 g, mag 16-bit, 100 Hz.
 * The sketch converts raw counts to the units the filter expects
 * (rad/s, m/s^2; magnetometer stays in raw-scaled uT).
 *
 * The magnetometer must be hard-iron calibrated for serious use:
 * subtract your offsets in readMag() where marked. Axis remap: the
 * AK8963 axes differ from the gyro/accel axes; the standard remap
 * (mx<->my swap, mz negated) is applied below so all three sensors
 * share one body frame.
 *
 * Output: roll, pitch, yaw in degrees plus trust flags at 20 Hz.
 */
#include <Wire.h>
#include <ArbiterAHRS.h>

static const uint8_t MPU = 0x68, MAG = 0x0C;
static const float DEG = 0.017453293f;
static const float GYRO_SCALE = 500.0f / 32768.0f * DEG;   /* rad/s per LSB */
static const float ACC_SCALE  = 4.0f * 9.80665f / 32768.0f; /* m/s^2 per LSB */

ArbiterAHRS ahrs;
float magAdj[3];                 /* AK8963 factory sensitivity adjustment */
uint32_t lastUs = 0, lastPrintMs = 0;

static void wr(uint8_t dev, uint8_t reg, uint8_t val) {
    Wire.beginTransmission(dev); Wire.write(reg); Wire.write(val);
    Wire.endTransmission();
}

static void rd(uint8_t dev, uint8_t reg, uint8_t n, uint8_t *buf) {
    Wire.beginTransmission(dev); Wire.write(reg);
    Wire.endTransmission(false);
    Wire.requestFrom(dev, n);
    for (uint8_t i = 0; i < n && Wire.available(); i++) buf[i] = Wire.read();
}

void setup() {
    Serial.begin(115200);
    Wire.begin();
    Wire.setClock(400000);

    wr(MPU, 0x6B, 0x01);         /* PWR_MGMT_1: wake, PLL clock */
    delay(50);
    wr(MPU, 0x1A, 0x03);         /* DLPF 41 Hz */
    wr(MPU, 0x1B, 0x08);         /* gyro +/-500 dps */
    wr(MPU, 0x1C, 0x08);         /* accel +/-4 g */
    wr(MPU, 0x1D, 0x03);         /* accel DLPF 41 Hz */
    wr(MPU, 0x37, 0x02);         /* INT_PIN_CFG: I2C bypass to AK8963 */

    uint8_t asa[3];
    wr(MAG, 0x0A, 0x00); delay(10);      /* power down */
    wr(MAG, 0x0A, 0x0F); delay(10);      /* fuse ROM access */
    rd(MAG, 0x10, 3, asa);
    for (int i = 0; i < 3; i++)
        magAdj[i] = ((float)asa[i] - 128.0f) / 256.0f + 1.0f;
    wr(MAG, 0x0A, 0x00); delay(10);
    wr(MAG, 0x0A, 0x16);                 /* 16-bit, continuous 100 Hz */

    ahrs.begin();
    lastUs = micros();
    Serial.println(F("MPU-9250 + ArbiterAHRS ready"));
}

static bool readMag(float m[3]) {
    uint8_t st, b[7];
    rd(MAG, 0x02, 1, &st);
    if (!(st & 0x01)) return false;      /* no new data */
    rd(MAG, 0x03, 7, b);                 /* HXL..HZH + ST2 */
    if (b[6] & 0x08) return false;       /* overflow */
    int16_t x = (int16_t)(b[1] << 8 | b[0]);
    int16_t y = (int16_t)(b[3] << 8 | b[2]);
    int16_t z = (int16_t)(b[5] << 8 | b[4]);
    /* remap AK8963 axes into the gyro/accel body frame */
    m[0] = (float)y * magAdj[1];         /* TODO: subtract hard-iron offset */
    m[1] = (float)x * magAdj[0];
    m[2] = -(float)z * magAdj[2];
    return true;
}

void loop() {
    uint8_t b[14];
    rd(MPU, 0x3B, 14, b);
    int16_t axr = (int16_t)(b[0] << 8 | b[1]);
    int16_t ayr = (int16_t)(b[2] << 8 | b[3]);
    int16_t azr = (int16_t)(b[4] << 8 | b[5]);
    int16_t gxr = (int16_t)(b[8] << 8 | b[9]);
    int16_t gyr = (int16_t)(b[10] << 8 | b[11]);
    int16_t gzr = (int16_t)(b[12] << 8 | b[13]);

    float mag[3];
    bool magOk = readMag(mag);

    uint32_t now = micros();
    float dt = (now - lastUs) * 1e-6f;
    lastUs = now;
    if (dt <= 0.0f || dt > 0.5f) return; /* first pass / timer wrap */

    if (magOk)
        ahrs.update(gxr*GYRO_SCALE, gyr*GYRO_SCALE, gzr*GYRO_SCALE,
                    axr*ACC_SCALE,  ayr*ACC_SCALE,  azr*ACC_SCALE,
                    mag[0], mag[1], mag[2], dt);
    else
        ahrs.updateIMU(gxr*GYRO_SCALE, gyr*GYRO_SCALE, gzr*GYRO_SCALE,
                       axr*ACC_SCALE,  ayr*ACC_SCALE,  azr*ACC_SCALE, dt);

    if (millis() - lastPrintMs >= 50) {
        lastPrintMs = millis();
        float r, p, y;
        ahrs.getEulerDeg(r, p, y);
        Serial.print(r, 1); Serial.print('\t');
        Serial.print(p, 1); Serial.print('\t');
        Serial.print(y, 1); Serial.print('\t');
        Serial.print(ahrs.accelTrusted() ? 'A' : '-');
        Serial.println(ahrs.magTrusted() ? 'M' : '-');
    }
}
