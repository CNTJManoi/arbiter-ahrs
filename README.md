# ARBITER firmware

ARBITER is an orientation filter (AHRS) for 9-axis MARG sensor arrays:
gyroscope, accelerometer, magnetometer. This directory contains the
platform-independent C99 core and ready-to-use packages for Arduino,
ESP-IDF, and STM32 HAL, each with examples.

## Algorithm in brief

The filter integrates the gyroscope and treats accelerometer and
magnetometer samples as requests for correction. A request passes three
consistency tests (specific-force norm, transport-theorem residual, and
magnetic norm/dip against online-learned references) and charges a
per-sensor evidence accumulator; corrections are enabled by a hysteresis
threshold on that accumulator. The magnitude of every correction is drawn
from a token bucket that refills at the expected gyroscope drift rate.
As a result, the total rotation that corrections can apply over any time
interval is bounded by the bucket capacity plus the integral of the refill
rate. A measurement that defeats all tests still cannot displace the
estimate faster than the gyroscope itself could drift.

Tilt corrections rotate about axes orthogonal to the vertical and heading
corrections rotate about the vertical, so a magnetometer fault cannot
corrupt roll and pitch. Gyro bias is measured directly during detected
rest, not accumulated from an error integral. The full derivation,
benchmark results (synthetic scenario and the 39-trial BROAD data set),
and the ablation study are in the accompanying paper draft in the parent
project.

Numbers that matter for firmware: 31 floating-point states plus eight binary flags (about 130 B), no
dynamic memory, no matrix operations, roughly 270 floating-point
operations plus 4 `sqrtf`, 3 `expf`, and 2 `atan2f` per step, float32
throughout.

## Directory layout

```
firmware/
├── core/                 platform-independent reference (arbiter.h/.c)
├── arduino/              Arduino library (library.properties format)
│   ├── src/              ArbiterAHRS.h wrapper + core copy
│   └── examples/
│       ├── NoHardware_Selftest/   runs without any sensor
│       └── Basic_MPU9250/         MPU-9250 over Wire (I2C)
├── esp32/                ESP-IDF component
│   ├── include/, arbiter.c, CMakeLists.txt, idf_component.yml
│   └── examples/basic/   self-test + optional MPU-9250 task
└── stm32/                core copy + HAL examples
    └── examples/
        ├── selftest_main.c        no sensor required, DWT timing
        └── mpu9250_f4_main.c      MPU-9250 on STM32F4 (CubeMX + HAL)
```

The four copies of `arbiter.h`/`arbiter.c` are identical; `core/` is the
source of truth. If you modify the core, copy it into the three platform
directories again.

## Input contract

| Input | Unit | Notes |
|---|---|---|
| gyroscope | rad/s | subtract nothing; the filter estimates bias itself |
| accelerometer | m/s² | specific force: +9.81 on Z when flat and still |
| magnetometer | any consistent unit | µT, gauss, or raw counts; the filter uses the direction and the relative norm |
| dt | s | measured interval since the previous call; jitter is acceptable |

Missing magnetometer samples are passed as `0,0,0` or `NAN`; the filter
then runs in 6-axis mode for that step and heading drifts at the residual
gyro rate. The magnetometer must be hard/soft-iron calibrated beforehand:
the filter rejects environmental disturbances, it does not calibrate the
sensor. Axes of all three sensors must form one right-handed body frame;
see the MPU-9250 examples for the AK8963 remap.

## Quick start

### Arduino (any core: AVR, SAMD, RP2040, ESP32-Arduino, STM32duino)

Copy `firmware/arduino/` into your sketchbook `libraries/` folder as
`ArbiterAHRS`, restart the IDE, and open
File → Examples → ArbiterAHRS → NoHardware_Selftest. The sketch prints
PASS when the filter converges on synthetic data; expected accuracy is
stated in the sketch header. `Basic_MPU9250` shows live fusion.

Minimal usage:

```cpp
#include <ArbiterAHRS.h>
ArbiterAHRS ahrs;

void setup() { ahrs.begin(); }

void loop() {
    // gx..gz rad/s, ax..az m/s^2, mx..mz any unit, dt seconds
    ahrs.update(gx, gy, gz, ax, ay, az, mx, my, mz, dt);
    float roll, pitch, yaw;
    ahrs.getEulerDeg(roll, pitch, yaw);
}
```

### ESP-IDF (ESP32, S2, S3, C3)

Copy `firmware/esp32/` into `<project>/components/arbiter` and add
`arbiter` to the `REQUIRES` of your main component, or build the bundled
example directly:

```sh
cd firmware/esp32/examples/basic
idf.py set-target esp32
idf.py build flash monitor
```

The example runs the self-test at boot and prints the measured per-step
time on your chip. Set `RUN_LIVE_MPU9250` to 1 in `main.c` to enable the
I2C task (SDA 21, SCL 22 by default).

### STM32 (HAL / CubeMX)

Add `stm32/arbiter.c` and `stm32/arbiter.h` to your CubeMX project. For a
first check, add `examples/selftest_main.c` and call `arbiter_selftest()`
after clock setup; it reports attitude, recovered bias, and the average
cycle count per step via the DWT counter. `examples/mpu9250_f4_main.c` is
a complete main loop for an MPU-9250 on I2C1 with UART output; replace
the marked lines with the init calls generated by CubeMX for your board.

## API (C core)

| Function | Purpose |
|---|---|
| `arbiter_default_params(ArbiterParams *p)` | fill `p` with the published defaults |
| `arbiter_init(ArbiterState *s, const ArbiterParams *p)` | reset state; `p` must outlive `s` |
| `arbiter_update(s, gx, gy, gz, ax, ay, az, mx, my, mz, dt)` | one fusion step |
| `arbiter_euler(s, &roll, &pitch, &yaw)` | ZYX Euler angles, rad |
| `s->q[4]` | quaternion, body to East-North-Up, Hamilton `[w x y z]` |
| `s->b[3]` | gyro bias estimate, rad/s |
| `s->gate_a`, `s->gate_m` | accelerometer / magnetometer currently trusted |
| `s->flag_dyn`, `s->flag_mag` | dynamic-acceleration / magnetic-disturbance flags |
| `s->is_static` | rest detector state |

The Arduino wrapper (`ArbiterAHRS.h`) exposes the same functionality as a
class with `update()`, `updateIMU()`, `getQuaternion()`, `getEulerDeg()`,
and diagnostic getters.

## Parameters

`arbiter_default_params()` sets values that performed well across the
synthetic benchmark and BROAD without per-dataset tuning. Two parameters
dominate accuracy:

| Parameter | Default | Meaning |
|---|---|---|
| `rho0` | 0.3 °/s | base refill rate of the correction budgets; raise it for gyros with worse bias stability |
| `diff_tol` | 0.8 rad/s | tolerance of the accelerometer transport-theorem test; lower it for cleaner mechanical setups |

The remaining fields are documented in `arbiter.h`. `cap_tilt` and
`cap_yaw` bound the single re-anchoring step after a distrust period;
`ev_up_tau_mag` sets how slowly magnetometer trust returns. One ablation
result is worth knowing: on data with mild disturbances, disabling the
hard gates (set `gate_hi` to 0 so gates stay open) and keeping only the
budget improved accuracy in our benchmark. Keep the gates for
environments with strong magnetic or inertial faults.

## Resource footprint

State: ~130 B RAM per filter instance plus ~100 B for the parameter
block. Code size: a few kB depending on the toolchain. Measured per-step
times are printed by the self-tests; the analytic estimates are 3–5 µs on
a Cortex-M7 at 216 MHz, 8–12 µs on a Cortex-M4F at 168 MHz, 15–25 µs on
an ESP32, 150–250 µs on a Cortex-M3 with software floats, and 1.5–2.5 ms
on an 8-bit AVR. On AVR, prefer the 6-axis mode and rates up to about
200 Hz.

## Validation

The self-test examples reproduce a static convergence case with known
ground truth: a fixed attitude, a 30-degree initial heading error, and an
injected gyro bias. The expected outcome is an attitude error below 1.5°
and a bias error below 0.15 °/s after 20 s at 250 Hz. Run the self-test
on your target before trusting live data: it exercises the full code path
and reports the real per-step cost on your hardware.

## Known limitations

Heading is unobservable without a magnetometer and drifts at the residual
bias rate in 6-axis mode. Gyro bias is measured only during detected
rest; applications that never rest accumulate the drift the budget is
sized for. Sustained coherent acceleration that keeps the specific-force
norm near g passes the diagnostics and steers tilt at the budgeted rate,
which is the designed worst case.


## Reproducibility package

The python/ directory contains the reference implementation and every script
behind the published results: the synthetic-scenario generator with fixed
seeds (simulation.py), the baseline filters (filters_baseline.py), the
reference ARBITER implementation (new_filter.py), the synthetic benchmark
with its tuning grids (benchmark.py), the BROAD benchmark under the official
TAGP protocol including VQF (broad_benchmark.py, vqf_wrapper.py), the
ablation study (ablation.py), figure generation (make_paper_figs.py), and
the result tables (python/results/). The uniform-sampling companion of the
synthetic scenario, which brings VQF into the stress test, is run by
benchmark_uniform.py (results in python/results/uniform_companion.json).
The per-trial BROAD errors behind the reported statistics are in
python/results/broad_per_trial.csv (39 trials by 6 filters; TAGP total,
heading, and inclination RMSE plus the zero-shot total); stats_broad.py
recomputes the Wilcoxon tests with Holm correction, the rank-biserial
effect sizes, and the Hodges-Lehmann estimates with confidence intervals.

To reproduce: install numpy, scipy, matplotlib, h5py, and vqf; clone the
BROAD data set (github.com/dlaidig/broad) into python/data/broad; then run
tests_sanity.py, benchmark.py, and broad_benchmark.py from the python/
directory. All algorithms are deterministic, so repeated runs reproduce the
tables exactly.
