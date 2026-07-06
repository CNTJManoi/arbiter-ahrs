/*
 * arbiter.h — ARBITER: Attitude by Rotation-Budgeted Independent
 * Triadic Evidence Reconciliation. Реализация для микроконтроллеров.
 *
 * Свойства:
 *   - float32, без динамической памяти, без внешних зависимостей (math.h);
 *   - состояние: 1 структура ~120 байт; стек update(): < 200 байт;
 *   - на шаг: ~260 FLOP + 4 sqrtf + 2-3 expf + 1-2 atan2f (см. arbiter.c);
 *   - пригодно для STM32 (HAL/LL), ESP32 (Arduino/IDF), Arduino Due;
 *     на AVR возможно на 50-100 Гц (программный float).
 *
 * Использование:
 *   ArbiterParams p;  ArbiterState s;
 *   arbiter_default_params(&p);
 *   arbiter_init(&s, &p);
 *   // в цикле измерений (СИ: рад/с, м/с^2, любые единицы магнитометра):
 *   arbiter_update(&s, gx, gy, gz, ax, ay, az, mx, my, mz, dt);
 *   // при отсутствии магнитометра передать mx=my=mz=0 (или NAN)
 *   float r, pch, yw;  arbiter_euler(&s, &r, &pch, &yw);
 *
 * Конвенции: земная СК ENU (x=восток, y=север, z=вверх); кватернион
 * q=[w,x,y,z] Hamilton, body->earth; углы Эйлера ZYX (roll,pitch,yaw).
 */
#ifndef ARBITER_H
#define ARBITER_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float g0;              /* модуль гравитации, м/с^2 (9.81) */
    /* токен-бакеты */
    float rho0;            /* базовая скорость пополнения бюджета, рад/с */
    float cap_tilt;        /* ёмкость бакета наклона, рад */
    float cap_yaw;         /* ёмкость бакета курса, рад */
    float routine_mult;    /* «дежурный» расход <= routine_mult*rho*dt */
    float burst_evidence;  /* burst из бакета только при E > этого порога */
    /* неопределённость дрейфа u_b */
    float ub0, ub_min, ub_max;   /* рад/с */
    float ub_grow;         /* рост, (рад/с)/с */
    float ub_decay_tau;    /* спад в статике, с */
    float ub_boost;        /* форсаж при насыщении, (рад/с)/с (эфф. x50) */
    float sat_thr_tilt, sat_thr_yaw;   /* пороги события насыщения, рад */
    /* накопители доказательств */
    float ev_up_tau, ev_up_tau_mag, ev_dn_tau;   /* с */
    float gate_hi, gate_lo;            /* гистерезис гейта */
    /* тесты акселерометра */
    float acc_norm_tol;    /* доля g */
    float diff_tol;        /* рад/с */
    float diff_leak_tau;   /* с */
    float omega_scale;     /* с: допуск *= (1 + omega_scale*|w|) */
    /* тесты магнитометра */
    float mag_norm_tol, dip_tol, mag_diff_tol, mag_hmin;
    /* статика и bias */
    float static_disp_thr; /* рад/с */
    float static_min_time; /* с */
    float lp_tau, bias_tau, bias_max;
    /* обучение эталонов поля */
    float field_tau;       /* с */
    float deadband;        /* рад */
} ArbiterParams;

typedef struct {
    const ArbiterParams *p;
    float q[4];            /* ориентация body->earth */
    float b[3];            /* оценка bias гироскопа, рад/с */
    float Ea, Em;          /* накопители доказательств [0..1] */
    bool  gate_a, gate_m;  /* состояния гейтов (гистерезис) */
    float Tt, Ty;          /* токен-бакеты, рад */
    float ub;              /* неопределённость дрейфа, рад/с */
    float B_ref, S_dip;    /* эталоны нормы поля и наклонения */
    bool  dip_init;
    float ca[3], cm[3];    /* накопители дифференциальных тестов */
    float pa[3], pm[3];    /* прошлые нормированные векторы */
    bool  pa_ok, pm_ok;
    float wl[3];           /* НЧ гироскопа */
    float disp;            /* НЧ |w - wl| */
    float static_t;
    bool  is_static;
    bool  initialized;
    /* диагностические флаги последнего шага */
    bool  flag_dyn;        /* динамическое ускорение (гейт acc закрыт) */
    bool  flag_mag;        /* магнитная помеха (гейт mag закрыт) */
} ArbiterState;

void arbiter_default_params(ArbiterParams *p);
void arbiter_init(ArbiterState *s, const ArbiterParams *p);
/* один шаг фильтра; при отсутствии магнитометра mx=my=mz=0 или NAN */
void arbiter_update(ArbiterState *s,
                    float gx, float gy, float gz,
                    float ax, float ay, float az,
                    float mx, float my, float mz,
                    float dt);
void arbiter_euler(const ArbiterState *s,
                   float *roll, float *pitch, float *yaw);

#ifdef __cplusplus
}
#endif
#endif /* ARBITER_H */
