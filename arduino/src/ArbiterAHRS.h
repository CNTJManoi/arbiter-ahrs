/*
 * ArbiterAHRS.h — Arduino wrapper for the ARBITER orientation filter.
 *
 * ARBITER fuses a 3-axis gyroscope, accelerometer, and (optionally)
 * magnetometer into a unit quaternion. Corrections from the vector sensors
 * are gated by per-sensor evidence accumulators and rate-limited by token
 * buckets whose refill rate tracks the expected gyroscope drift, so a
 * corrupted measurement cannot displace the estimate faster than the
 * gyroscope itself could drift.
 *
 * Input units (fixed contract):
 *   gyroscope      rad/s
 *   accelerometer  m/s^2  (specific force; +9.81 on Z when flat and still)
 *   magnetometer   any consistent unit (uT, gauss, raw counts);
 *                  pass 0,0,0 or NAN when the sample is unavailable
 *   dt             seconds since the previous update
 *
 * The magnetometer must be hard/soft-iron calibrated beforehand. The filter
 * handles environmental disturbances; it does not calibrate the sensor.
 */
#ifndef ARBITER_AHRS_ARDUINO_H
#define ARBITER_AHRS_ARDUINO_H

#include <Arduino.h>
extern "C" {
#include "arbiter.h"
}

class ArbiterAHRS {
public:
    ArbiterAHRS() { arbiter_default_params(&params_); }

    /* Call once, after adjusting parameters if needed. */
    void begin() { arbiter_init(&state_, &params_); }

    /* Reset the filter state; parameters are kept. */
    void reset() { arbiter_init(&state_, &params_); }

    /* One fusion step. Returns nothing; read the results via getters. */
    void update(float gx, float gy, float gz,
                float ax, float ay, float az,
                float mx, float my, float mz,
                float dt) {
        arbiter_update(&state_, gx, gy, gz, ax, ay, az, mx, my, mz, dt);
    }

    /* 6-axis variant (no magnetometer): heading is stabilized by the
       gyroscope only and will drift slowly. */
    void updateIMU(float gx, float gy, float gz,
                   float ax, float ay, float az, float dt) {
        arbiter_update(&state_, gx, gy, gz, ax, ay, az,
                       0.0f, 0.0f, 0.0f, dt);
    }

    /* Orientation, body frame to East-North-Up, Hamilton convention. */
    void getQuaternion(float &w, float &x, float &y, float &z) const {
        w = state_.q[0]; x = state_.q[1]; y = state_.q[2]; z = state_.q[3];
    }

    /* Euler angles in radians (roll X, pitch Y, yaw Z, ZYX sequence). */
    void getEuler(float &roll, float &pitch, float &yaw) const {
        arbiter_euler(&state_, &roll, &pitch, &yaw);
    }

    /* Euler angles in degrees. */
    void getEulerDeg(float &roll, float &pitch, float &yaw) const {
        arbiter_euler(&state_, &roll, &pitch, &yaw);
        roll *= 57.29578f; pitch *= 57.29578f; yaw *= 57.29578f;
    }

    /* Diagnostics. */
    bool  accelTrusted()  const { return state_.gate_a; }
    bool  magTrusted()    const { return state_.gate_m; }
    bool  isStatic()      const { return state_.is_static; }
    bool  dynAccelFlag()  const { return state_.flag_dyn; }
    bool  magDisturbFlag()const { return state_.flag_mag; }
    float accelEvidence() const { return state_.Ea; }
    float magEvidence()   const { return state_.Em; }
    void  getGyroBias(float &bx, float &by, float &bz) const {
        bx = state_.b[0]; by = state_.b[1]; bz = state_.b[2];
    }

    /* Direct access for advanced use (parameter tuning before begin()). */
    ArbiterParams &parameters() { return params_; }
    const ArbiterState &state() const { return state_; }

private:
    ArbiterParams params_;
    ArbiterState  state_;
};

#endif /* ARBITER_AHRS_ARDUINO_H */
