/*
 * arbiter.c — ARBITER для микроконтроллеров (float32, C99).
 * Алгоритм и параметры соответствуют описанию в сопроводительной статье.
 *
 * Бюджет вычислений на шаг (примерно, с магнитометром):
 *   ~270 FLOP (умножения/сложения) + 4 sqrtf + 3 expf + 2 atan2f + 2 sincosf.
 *   На Cortex-M4F (STM32F407, 168 МГц, FPU): ~8-12 мкс  -> >80 кГц запас;
 *   на Cortex-M7  (STM32F767): ~3-5 мкс; на ESP32: ~15-25 мкс;
 *   на Cortex-M3 без FPU (STM32F103, программный float): ~150-250 мкс;
 *   на AVR (16 МГц): ~1.5-2.5 мс -> максимум ~100-200 Гц без магнитометра.
 */
#include "arbiter.h"
#include <math.h>

#define DEG2RAD 0.017453292519943295f

/* ------------------------------------------------- кватернионные операции */
static void quat_mult(const float a[4], const float b[4], float out[4])
{
    out[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    out[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
    out[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
    out[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
}

static void quat_normalize(float q[4])
{
    float n = sqrtf(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    if (n < 1e-12f) { q[0] = 1.0f; q[1] = q[2] = q[3] = 0.0f; return; }
    float inv = 1.0f / n;
    q[0] *= inv; q[1] *= inv; q[2] *= inv; q[3] *= inv;
}

static void quat_from_rotvec(float rx, float ry, float rz, float out[4])
{
    float ang = sqrtf(rx*rx + ry*ry + rz*rz);
    if (ang < 1e-9f) {
        out[0] = 1.0f; out[1] = 0.5f*rx; out[2] = 0.5f*ry; out[3] = 0.5f*rz;
        quat_normalize(out);
        return;
    }
    float half = 0.5f * ang;
    float s = sinf(half) / ang;
    out[0] = cosf(half); out[1] = rx*s; out[2] = ry*s; out[3] = rz*s;
}

/* v_e = R(q) v_b */
static void quat_rotate(const float q[4], const float v[3], float out[3])
{
    float tx = 2.0f*(q[2]*v[2] - q[3]*v[1]);
    float ty = 2.0f*(q[3]*v[0] - q[1]*v[2]);
    float tz = 2.0f*(q[1]*v[1] - q[2]*v[0]);
    out[0] = v[0] + q[0]*tx + (q[2]*tz - q[3]*ty);
    out[1] = v[1] + q[0]*ty + (q[3]*tx - q[1]*tz);
    out[2] = v[2] + q[0]*tz + (q[1]*ty - q[2]*tx);
}

/* R(q)^T e_z — ожидаемое направление гравитации в body */
static void gravity_body(const float q[4], float out[3])
{
    out[0] = 2.0f*(q[1]*q[3] - q[0]*q[2]);
    out[1] = 2.0f*(q[2]*q[3] + q[0]*q[1]);
    out[2] = 1.0f - 2.0f*(q[1]*q[1] + q[2]*q[2]);
}

static float gauss_score(float x, float tol)
{
    float r = x / tol;
    return expf(-r * r);
}

/* --------------------------------------------------------------- params */
void arbiter_default_params(ArbiterParams *p)
{
    p->g0 = 9.81f;
    p->rho0 = 0.3f * DEG2RAD;
    p->cap_tilt = 10.0f * DEG2RAD;
    p->cap_yaw = 15.0f * DEG2RAD;
    p->routine_mult = 2.0f;
    p->burst_evidence = -1.0f;   /* двухскоростная политика отключена */
    p->ub0 = 2.0f * DEG2RAD;
    p->ub_min = 0.05f * DEG2RAD;
    p->ub_max = 2.0f * DEG2RAD;
    p->ub_grow = 0.005f * DEG2RAD;
    p->ub_decay_tau = 3.0f;
    p->ub_boost = 0.2f * DEG2RAD;
    p->sat_thr_tilt = 1.5f * DEG2RAD;
    p->sat_thr_yaw = 3.0f * DEG2RAD;
    p->ev_up_tau = 0.3f;
    p->ev_up_tau_mag = 2.0f;   /* курс: медленное возвращение доверия */
    p->ev_dn_tau = 0.04f;
    p->gate_hi = 0.35f;
    p->gate_lo = 0.20f;
    p->acc_norm_tol = 0.12f;
    p->diff_tol = 0.8f;
    p->diff_leak_tau = 0.15f;
    p->omega_scale = 0.5f;
    p->mag_norm_tol = 0.15f;
    p->dip_tol = 0.12f;
    p->mag_diff_tol = 0.8f;
    p->mag_hmin = 0.25f;
    p->static_disp_thr = 0.3f * DEG2RAD;
    p->static_min_time = 0.5f;
    p->lp_tau = 0.5f;
    p->bias_tau = 2.0f;
    p->bias_max = 3.0f * DEG2RAD;
    p->field_tau = 30.0f;
    p->deadband = 0.02f * DEG2RAD;
}

void arbiter_init(ArbiterState *s, const ArbiterParams *p)
{
    int i;
    s->p = p;
    s->q[0] = 1.0f; s->q[1] = s->q[2] = s->q[3] = 0.0f;
    for (i = 0; i < 3; i++) {
        s->b[i] = 0.0f; s->ca[i] = 0.0f; s->cm[i] = 0.0f;
        s->pa[i] = 0.0f; s->pm[i] = 0.0f; s->wl[i] = 0.0f;
    }
    s->Ea = 0.0f; s->Em = 0.0f;
    s->gate_a = false; s->gate_m = false;
    s->Tt = p->cap_tilt; s->Ty = p->cap_yaw;   /* бакеты стартуют полными */
    s->ub = p->ub0;
    s->B_ref = -1.0f; s->S_dip = 0.0f; s->dip_init = false;
    s->pa_ok = false; s->pm_ok = false;
    s->disp = 1.0f; s->static_t = 0.0f; s->is_static = false;
    s->initialized = false;
    s->flag_dyn = false; s->flag_mag = false;
}

/* инициализация ориентации по первому измерению acc (+mag) */
static void first_fix(ArbiterState *s,
                      float ax, float ay, float az,
                      float mx, float my, float mz, bool has_mag)
{
    float na = sqrtf(ax*ax + ay*ay + az*az);
    if (na < 1e-6f) return;
    ax /= na; ay /= na; az /= na;
    /* тилт: повернуть a-hat в e_z; ось = a x e_z = (ay, -ax, 0) */
    float sn = sqrtf(ax*ax + ay*ay);
    if (sn < 1e-9f) {
        s->q[0] = (az > 0.0f) ? 1.0f : 0.0f;
        s->q[1] = (az > 0.0f) ? 0.0f : 1.0f;
        s->q[2] = s->q[3] = 0.0f;
    } else {
        float ang = atan2f(sn, az);
        quat_from_rotvec(ay/sn*ang, -ax/sn*ang, 0.0f, s->q);
    }
    if (has_mag) {
        float m[3] = { mx, my, mz }, me[3];
        quat_rotate(s->q, m, me);
        float h = sqrtf(me[0]*me[0] + me[1]*me[1]);
        if (h > 1e-6f) {
            float psi = atan2f(me[0], me[1]);   /* ошибка yaw = -psi */
            float qz[4], tmp[4];
            quat_from_rotvec(0.0f, 0.0f, psi, qz);
            quat_mult(qz, s->q, tmp);
            tmp[0] = tmp[0]; /* copy */
            s->q[0] = tmp[0]; s->q[1] = tmp[1];
            s->q[2] = tmp[2]; s->q[3] = tmp[3];
        }
        s->B_ref = sqrtf(mx*mx + my*my + mz*mz);
        if (s->B_ref > 1e-6f) {
            s->pm[0] = mx / s->B_ref; s->pm[1] = my / s->B_ref;
            s->pm[2] = mz / s->B_ref;
            s->pm_ok = true;
        }
    }
    s->pa[0] = ax; s->pa[1] = ay; s->pa[2] = az; s->pa_ok = true;
    quat_normalize(s->q);
    s->initialized = true;
}

void arbiter_update(ArbiterState *s,
                    float gx, float gy, float gz,
                    float ax, float ay, float az,
                    float mx, float my, float mz, float dt)
{
    const ArbiterParams *p = s->p;
    bool has_mag = (mx == mx) && (my == my) && (mz == mz) &&
                   (mx != 0.0f || my != 0.0f || mz != 0.0f);
    if (!s->initialized) {
        first_fix(s, ax, ay, az, mx, my, mz, has_mag);
        return;
    }

    /* ===== 1. Предсказание гироскопом */
    float wx = gx - s->b[0], wy = gy - s->b[1], wz = gz - s->b[2];
    float dq[4], qn[4];
    quat_from_rotvec(wx*dt, wy*dt, wz*dt, dq);
    quat_mult(s->q, dq, qn);

    /* ===== 2. Пополнение бюджетов */
    float rho = p->rho0 + s->ub;
    s->Tt += rho * dt; if (s->Tt > p->cap_tilt) s->Tt = p->cap_tilt;
    s->Ty += rho * dt; if (s->Ty > p->cap_yaw) s->Ty = p->cap_yaw;
    s->ub += p->ub_grow * dt; if (s->ub > p->ub_max) s->ub = p->ub_max;

    float wnorm = sqrtf(wx*wx + wy*wy + wz*wz);
    float tol_scale = 1.0f + p->omega_scale * wnorm;

    /* ===== 3. Диагностика акселерометра */
    float na = sqrtf(ax*ax + ay*ay + az*az);
    float s_a = 0.0f;
    float ah[3] = { 0.0f, 0.0f, 1.0f };
    if (na > 1e-6f) {
        ah[0] = ax/na; ah[1] = ay/na; ah[2] = az/na;
        float s_norm = gauss_score(na/p->g0 - 1.0f, p->acc_norm_tol);
        float lam = expf(-dt / p->diff_leak_tau);
        if (s->pa_ok) {
            s->ca[0] = s->ca[0]*lam + (ah[0]-s->pa[0]) + (wy*ah[2]-wz*ah[1])*dt;
            s->ca[1] = s->ca[1]*lam + (ah[1]-s->pa[1]) + (wz*ah[0]-wx*ah[2])*dt;
            s->ca[2] = s->ca[2]*lam + (ah[2]-s->pa[2]) + (wx*ah[1]-wy*ah[0])*dt;
        }
        float cn = sqrtf(s->ca[0]*s->ca[0] + s->ca[1]*s->ca[1]
                         + s->ca[2]*s->ca[2]) / p->diff_leak_tau;
        float s_diff = gauss_score(cn, p->diff_tol * tol_scale);
        s_a = s_norm * s_diff;
        s->pa[0] = ah[0]; s->pa[1] = ah[1]; s->pa[2] = ah[2]; s->pa_ok = true;
    }
    {   /* накопитель доказательств (асимметричный) */
        float tau = (s_a > s->Ea) ? p->ev_up_tau : p->ev_dn_tau;
        s->Ea += (s_a - s->Ea) * (dt / (tau + dt));
        if (s->Ea < 0.0f) s->Ea = 0.0f;
        if (s->Ea > 1.0f) s->Ea = 1.0f;
    }
    s->gate_a = s->gate_a ? (s->Ea > p->gate_lo) : (s->Ea > p->gate_hi);
    s->flag_dyn = !s->gate_a;

    /* ===== 4. Квазистатика и оценка bias */
    {
        float k_lp = dt / (p->lp_tau + dt);
        s->wl[0] += (gx - s->wl[0]) * k_lp;
        s->wl[1] += (gy - s->wl[1]) * k_lp;
        s->wl[2] += (gz - s->wl[2]) * k_lp;
        float dx = gx - s->wl[0], dy = gy - s->wl[1], dz = gz - s->wl[2];
        float dev = sqrtf(dx*dx + dy*dy + dz*dz);
        s->disp += (dev - s->disp) * k_lp;
        bool quasi = (s->disp < p->static_disp_thr) && (s_a > 0.5f);
        s->static_t = quasi ? (s->static_t + dt) : 0.0f;
        s->is_static = (s->static_t >= p->static_min_time);
        if (s->is_static) {
            float k_b = dt / (p->bias_tau + dt);
            int i;
            for (i = 0; i < 3; i++) {
                s->b[i] += (s->wl[i] - s->b[i]) * k_b;
                if (s->b[i] > p->bias_max) s->b[i] = p->bias_max;
                if (s->b[i] < -p->bias_max) s->b[i] = -p->bias_max;
            }
            s->ub *= expf(-dt / p->ub_decay_tau);
            if (s->ub < p->ub_min) s->ub = p->ub_min;
        }
    }

    /* ===== 5. Коррекция наклона: геодезический снап в пределах бюджета */
    bool sat_tilt = false;
    if (s->gate_a && na > 1e-6f) {
        float v[3];
        gravity_body(qn, v);
        float ex = ah[1]*v[2] - ah[2]*v[1];
        float ey = ah[2]*v[0] - ah[0]*v[2];
        float ez = ah[0]*v[1] - ah[1]*v[0];
        float sn = sqrtf(ex*ex + ey*ey + ez*ez);
        float cs = ah[0]*v[0] + ah[1]*v[1] + ah[2]*v[2];
        float theta = atan2f(sn, cs);
        if (theta > p->deadband && sn > 1e-9f) {
            float allow = s->Tt;
            if (s->Ea < p->burst_evidence) {   /* двухскоростная политика */
                float routine = p->routine_mult * rho * dt;
                if (routine < allow) allow = routine;
            }
            float theta_c = (theta < allow) ? theta : allow;
            s->Tt -= theta_c;
            float f = theta_c / sn;
            float qc[4], tmp[4];
            quat_from_rotvec(ex*f, ey*f, ez*f, qc);
            quat_mult(qn, qc, tmp);
            qn[0]=tmp[0]; qn[1]=tmp[1]; qn[2]=tmp[2]; qn[3]=tmp[3];
            if (theta - theta_c > p->sat_thr_tilt) sat_tilt = true;
        }
    }

    /* ===== 6. Диагностика магнитометра */
    float s_m = 0.0f;
    float mh[3] = { 0.0f, 0.0f, 0.0f };
    float nm = 0.0f;
    if (has_mag) {
        nm = sqrtf(mx*mx + my*my + mz*mz);
        if (nm > 1e-6f) {
            mh[0] = mx/nm; mh[1] = my/nm; mh[2] = mz/nm;
            if (s->B_ref <= 0.0f) s->B_ref = nm;
            float s_norm_m = gauss_score(nm/s->B_ref - 1.0f, p->mag_norm_tol);
            float lam = expf(-dt / p->diff_leak_tau);
            if (s->pm_ok) {
                s->cm[0] = s->cm[0]*lam + (mh[0]-s->pm[0]) + (wy*mh[2]-wz*mh[1])*dt;
                s->cm[1] = s->cm[1]*lam + (mh[1]-s->pm[1]) + (wz*mh[0]-wx*mh[2])*dt;
                s->cm[2] = s->cm[2]*lam + (mh[2]-s->pm[2]) + (wx*mh[1]-wy*mh[0])*dt;
            }
            float cnm = sqrtf(s->cm[0]*s->cm[0] + s->cm[1]*s->cm[1]
                              + s->cm[2]*s->cm[2]) / p->diff_leak_tau;
            float s_diff_m = gauss_score(cnm, p->mag_diff_tol * tol_scale);
            float dv[3];
            if (na > 1e-6f && s->gate_a) {
                dv[0] = ah[0]; dv[1] = ah[1]; dv[2] = ah[2];
            } else {
                gravity_body(qn, dv);
            }
            float dip = dv[0]*mh[0] + dv[1]*mh[1] + dv[2]*mh[2];
            if (!s->dip_init) { s->S_dip = dip; s->dip_init = true; }
            float s_dip = gauss_score(dip - s->S_dip, p->dip_tol);
            s_m = s_norm_m * s_diff_m * s_dip;
            s->pm[0]=mh[0]; s->pm[1]=mh[1]; s->pm[2]=mh[2]; s->pm_ok = true;
            if (s_m > 0.7f && s->gate_a) {   /* онлайн-обучение эталонов */
                float k_f = dt / (p->field_tau + dt);
                s->B_ref += (nm - s->B_ref) * k_f;
                s->S_dip += (dip - s->S_dip) * k_f;
            }
        }
    } else {
        s->pm_ok = false;
        s->cm[0] = s->cm[1] = s->cm[2] = 0.0f;
    }
    {
        float tau = (s_m > s->Em) ? p->ev_up_tau_mag : p->ev_dn_tau;
        s->Em += (s_m - s->Em) * (dt / (tau + dt));
        if (s->Em < 0.0f) s->Em = 0.0f;
        if (s->Em > 1.0f) s->Em = 1.0f;
    }
    s->gate_m = s->gate_m ? (s->Em > p->gate_lo) : (s->Em > p->gate_hi);
    s->flag_mag = has_mag && !s->gate_m;

    /* ===== 7. Коррекция yaw вокруг земной вертикали в пределах бюджета */
    bool sat_yaw = false;
    if (has_mag && s->gate_m && s->gate_a && nm > 1e-6f) {
        float me[3];
        quat_rotate(qn, mh, me);
        float h = sqrtf(me[0]*me[0] + me[1]*me[1]);
        if (h > p->mag_hmin) {
            float psi = atan2f(me[0], me[1]);   /* ошибка yaw = -psi */
            float apsi = (psi >= 0.0f) ? psi : -psi;
            if (apsi > p->deadband) {
                float allow = s->Ty;
                if (s->Em < p->burst_evidence) {
                    float routine = p->routine_mult * rho * dt;
                    if (routine < allow) allow = routine;
                }
                float psi_c = (apsi < allow) ? apsi : allow;
                s->Ty -= psi_c;
                float sgn = (psi >= 0.0f) ? 1.0f : -1.0f;
                float qz[4], tmp[4];
                quat_from_rotvec(0.0f, 0.0f, sgn*psi_c, qz);
                quat_mult(qz, qn, tmp);
                qn[0]=tmp[0]; qn[1]=tmp[1]; qn[2]=tmp[2]; qn[3]=tmp[3];
                if (apsi - psi_c > p->sat_thr_yaw) sat_yaw = true;
            }
        }
    }

    /* ===== 8. Насыщение бюджета -> форсаж неопределённости */
    if (sat_tilt || sat_yaw) {
        s->ub += p->ub_boost * dt * 10.0f;
        if (s->ub > p->ub_max) s->ub = p->ub_max;
    }

    quat_normalize(qn);
    s->q[0]=qn[0]; s->q[1]=qn[1]; s->q[2]=qn[2]; s->q[3]=qn[3];
}

void arbiter_euler(const ArbiterState *s,
                   float *roll, float *pitch, float *yaw)
{
    const float *q = s->q;
    float sinp = 2.0f*(q[0]*q[2] - q[3]*q[1]);
    if (sinp > 1.0f) sinp = 1.0f;
    if (sinp < -1.0f) sinp = -1.0f;
    *roll = atan2f(2.0f*(q[0]*q[1] + q[2]*q[3]),
                   1.0f - 2.0f*(q[1]*q[1] + q[2]*q[2]));
    *pitch = asinf(sinp);
    *yaw = atan2f(2.0f*(q[0]*q[3] + q[1]*q[2]),
                  1.0f - 2.0f*(q[2]*q[2] + q[3]*q[3]));
}
