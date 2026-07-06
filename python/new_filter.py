"""
new_filter.py — ARBITER: Attitude by Rotation-Budgeted Independent
Triadic Evidence Reconciliation. Кандидат нового алгоритма ориентации.

Три несущих механизма (отсутствуют в Madgwick/Mahony/Complementary/EKF):

1. ТОКЕН-БАКЕТ БЮДЖЕТА КОРРЕКЦИИ. Любая коррекция кватерниона тратит
   скалярный «бюджет поворота» (радианы). Бюджет пополняется со скоростью,
   равной физически ожидаемому дрейфу гироскопа rho = rho0 + u_b, где u_b —
   скалярная неопределённость bias (растёт со временем, падает при статике,
   форсируется при устойчивом насыщении бюджета). Ёмкость ограничена, поэтому
   после длительного недоверия допускается ограниченный «снап», а любой
   сенсорный обман не может увести оценку быстрее, чем гироскоп мог бы
   реально задрейфовать. Раздельные бакеты для наклона (T_t) и yaw (T_y).

2. АСИММЕТРИЧНЫЕ НАКОПИТЕЛИ ДОКАЗАТЕЛЬСТВ E_a, E_m ∈ [0,1] с гистерезисным
   жёстким гейтом: согласованные измерения заряжают медленно (tau_up),
   нарушения разряжают быстро (tau_dn << tau_up). Коррекция ВКЛ/ВЫКЛ
   (событийная), а не непрерывно взвешенная.

3. ДИФФЕРЕНЦИАЛЬНЫЕ ТЕСТЫ СОГЛАСОВАННОСТИ. Для единичного вектора v̂,
   измеренного в body-СК, при чистом вращении dv̂/dt = −ω × v̂. Накопитель
   c ← λ·c + (Δv̂ + ω̂×v̂·dt) отделяет: для акселерометра — динамическое
   ускорение, для магнитометра — изменение самого поля (помеху) независимо
   от вращения. Плюс тесты нормы поля и магнитного наклонения (dip),
   эталоны которых обучаются онлайн в доверенные периоды.

Коррекции — дискретные геодезические шаги «снап в пределах бюджета»
theta_c = min(theta_err, T): НЕ пропорциональная обратная связь (Mahony),
НЕ градиентный шаг (Madgwick), НЕ ковариационное взвешивание (EKF).
Bias гироскопа оценивается ТОЛЬКО в валидированные квазистатические эпизоды
прямым усреднением показаний гироскопа (без интегрального члена Mahony).

Память O(1), матриц нет, обращений нет, один sqrt/atan2 на канал.
"""
from math import sqrt, atan2, exp
from imu_math import (quat_mult, quat_normalize, quat_from_rotvec,
                      quat_rotate, gravity_body)
from filters_baseline import init_quat_from_accel_mag, _mag_ok

DEG = 3.14159265358979 / 180.0


class ArbiterFilter:
    name = 'ARBITER'

    def __init__(self,
                 g0=9.81,
                 # --- токен-бакеты (двухскоростная политика расхода:
                 # «дежурный» расход <= routine_mult*rho*dt всегда при
                 # открытых воротах; burst из накопленного бакета — только
                 # при высоком доверии E > burst_evidence)
                 rho0=0.3*DEG,            # базовая скорость пополнения, рад/с
                 cap_tilt=10.0*DEG, cap_yaw=15.0*DEG,
                 # burst_evidence=-1: двухскоростная политика ОТКЛЮЧЕНА.
                 # Проверена и отвергнута: душит перепривязку после потери
                 # магнитометра (см. results/ablation.md, вариант dualrate)
                 routine_mult=2.0, burst_evidence=-1.0,
                 # --- неопределённость дрейфа u_b (управляет пополнением)
                 ub0=2.0*DEG, ub_min=0.05*DEG, ub_max=2.0*DEG,
                 ub_grow=0.005*DEG,       # рост, (рад/с)/с
                 ub_decay_tau=3.0,        # спад при статике, с
                 ub_boost=0.2*DEG,        # форсаж при насыщении, (рад/с)/с
                 sat_thr_tilt=1.5*DEG, sat_thr_yaw=3.0*DEG,
                 # --- накопители доказательств
                 # доверие акселерометру возвращается быстро (наклон
                 # восстановим), магнитометру — медленно (курс по гравитации
                 # невосстановим, ложное переоткрытие ворот опасно)
                 ev_up_tau=0.3, ev_up_tau_mag=2.0, ev_dn_tau=0.04,
                 gate_hi=0.35, gate_lo=0.20,
                 # --- тесты акселерометра
                 acc_norm_tol=0.12,       # доля g
                 diff_tol=0.8,            # рад/с (базовый допуск)
                 diff_leak_tau=0.15,      # с, утечка накопителя c
                 omega_scale=0.5,         # с: допуск x(1+omega_scale*|w|)
                 # --- тесты магнитометра
                 mag_norm_tol=0.15, dip_tol=0.12, mag_diff_tol=0.8,
                 mag_hmin=0.25,
                 # --- статика и bias
                 static_disp_thr=0.3*DEG, static_min_time=0.5,
                 lp_tau=0.5, bias_tau=2.0, bias_max=3.0*DEG,
                 # --- обучение эталонов поля
                 field_tau=30.0,
                 deadband=0.02*DEG):
        p = dict(locals()); p.pop('self')
        self.p = p
        self.reset()

    def reset(self):
        p = self.p
        self.q = (1.0, 0.0, 0.0, 0.0)
        self.bx = self.by = self.bz = 0.0
        self.Ea = 0.0
        self.Em = 0.0
        self.gate_a = False
        self.gate_m = False
        self.Tt = p['cap_tilt']          # бакеты стартуют полными:
        self.Ty = p['cap_yaw']           # быстрая начальная сходимость
        self.ub = p['ub0']
        self.B_ref = -1.0                # эталон нормы поля (обучается)
        self.S_dip = 0.0                 # эталон â·m̂ = −sin(dip) (обучается)
        self._dip_init = False
        # накопители дифференциальных тестов
        self.cax = self.cay = self.caz = 0.0
        self.cmx = self.cmy = self.cmz = 0.0
        self._pa = None                  # прошлое â
        self._pm = None                  # прошлое m̂
        # детектор статики
        self.wlx = self.wly = self.wlz = 0.0   # НЧ гироскопа
        self.disp = 1.0                  # НЧ |ω − ω_lp|
        self.static_t = 0.0
        self.is_static = False
        self._init = False

    # ------------------------------------------------------------------ util
    @staticmethod
    def _gauss(x, tol):
        r = x / tol
        return exp(-r*r)

    def _evidence(self, E, s, dt, up_tau=None):
        p = self.p
        tau = (up_tau if up_tau is not None else p['ev_up_tau']) \
            if s > E else p['ev_dn_tau']
        E += (s - E) * (dt / (tau + dt))
        return 0.0 if E < 0.0 else (1.0 if E > 1.0 else E)

    def _gate(self, gate, E):
        p = self.p
        if gate:
            return E > p['gate_lo']
        return E > p['gate_hi']

    # ---------------------------------------------------------------- update
    def update(self, gyr, acc, mag, dt):
        p = self.p
        has_mag = _mag_ok(mag)
        if not self._init:
            self.q = init_quat_from_accel_mag(acc, mag if has_mag else None)
            na = sqrt(acc[0]**2 + acc[1]**2 + acc[2]**2)
            if na > 1e-6:
                self._pa = (acc[0]/na, acc[1]/na, acc[2]/na)
            if has_mag:
                nm = sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
                self.B_ref = nm
                if nm > 1e-6:
                    self._pm = (mag[0]/nm, mag[1]/nm, mag[2]/nm)
            self._init = True
            return self.q

        # ===== 1. Гироскопический наблюдатель: предсказание
        wx = gyr[0] - self.bx
        wy = gyr[1] - self.by
        wz = gyr[2] - self.bz
        q = quat_mult(self.q, quat_from_rotvec(wx*dt, wy*dt, wz*dt))

        # ===== 2. Пополнение бюджетов (rho = rho0 + u_b)
        rho = p['rho0'] + self.ub
        self.Tt = min(self.Tt + rho*dt, p['cap_tilt'])
        self.Ty = min(self.Ty + rho*dt, p['cap_yaw'])
        self.ub = min(self.ub + p['ub_grow']*dt, p['ub_max'])

        # ===== 3. Диагностика акселерометра
        ax, ay, az = acc
        na = sqrt(ax*ax + ay*ay + az*az)
        wnorm = sqrt(wx*wx + wy*wy + wz*wz)
        tol_scale = 1.0 + p['omega_scale'] * wnorm   # допуски растут с |w|:
        s_a = 0.0                                    # дискретизация/рассинхрон
        dyn_flag = True
        if na > 1e-6:
            ah = (ax/na, ay/na, az/na)
            # тест нормы
            s_norm = self._gauss(na/p['g0'] - 1.0, p['acc_norm_tol'])
            # дифференциальный тест: c ← λc + (Δâ + ω̂×â·dt)
            lam = exp(-dt / p['diff_leak_tau'])
            if self._pa is not None:
                dx_ = ah[0] - self._pa[0] + (wy*ah[2] - wz*ah[1])*dt
                dy_ = ah[1] - self._pa[1] + (wz*ah[0] - wx*ah[2])*dt
                dz_ = ah[2] - self._pa[2] + (wx*ah[1] - wy*ah[0])*dt
                self.cax = self.cax*lam + dx_
                self.cay = self.cay*lam + dy_
                self.caz = self.caz*lam + dz_
            cn = sqrt(self.cax**2 + self.cay**2 + self.caz**2) / p['diff_leak_tau']
            s_diff = self._gauss(cn, p['diff_tol'] * tol_scale)
            s_a = s_norm * s_diff
            self._pa = ah
        self.Ea = self._evidence(self.Ea, s_a, dt)
        self.gate_a = self._gate(self.gate_a, self.Ea)
        dyn_flag = not self.gate_a

        # ===== 4. Детектор квазистатики и оценка bias (только в статике)
        k_lp = dt / (p['lp_tau'] + dt)
        self.wlx += (gyr[0] - self.wlx) * k_lp
        self.wly += (gyr[1] - self.wly) * k_lp
        self.wlz += (gyr[2] - self.wlz) * k_lp
        dev = sqrt((gyr[0]-self.wlx)**2 + (gyr[1]-self.wly)**2
                   + (gyr[2]-self.wlz)**2)
        self.disp += (dev - self.disp) * k_lp
        quasi = (self.disp < p['static_disp_thr']) and (s_a > 0.5)
        if quasi:
            self.static_t += dt
        else:
            self.static_t = 0.0
        self.is_static = self.static_t >= p['static_min_time']
        if self.is_static:
            k_b = dt / (p['bias_tau'] + dt)
            bm = p['bias_max']
            self.bx = max(-bm, min(bm, self.bx + (self.wlx - self.bx)*k_b))
            self.by = max(-bm, min(bm, self.by + (self.wly - self.by)*k_b))
            self.bz = max(-bm, min(bm, self.bz + (self.wlz - self.bz)*k_b))
            self.ub = max(p['ub_min'], self.ub * exp(-dt / p['ub_decay_tau']))

        # ===== 5. Коррекция наклона: геодезический снап в пределах бюджета
        sat_tilt = False
        if self.gate_a and na > 1e-6:
            vx, vy, vz = gravity_body(q)
            ex = ah[1]*vz - ah[2]*vy
            ey = ah[2]*vx - ah[0]*vz
            ez = ah[0]*vy - ah[1]*vx          # â × v̂
            sn = sqrt(ex*ex + ey*ey + ez*ez)
            cs = ah[0]*vx + ah[1]*vy + ah[2]*vz
            theta = atan2(sn, cs)
            if theta > p['deadband'] and sn > 1e-9:
                allow = self.Tt
                if self.Ea < p['burst_evidence']:
                    # пограничное доверие: только «дежурный» расход,
                    # burst из накопленного бакета не разрешён
                    routine = p['routine_mult'] * rho * dt
                    allow = allow if allow < routine else routine
                theta_c = theta if theta < allow else allow
                self.Tt -= theta_c
                f = theta_c / sn
                q = quat_mult(q, quat_from_rotvec(ex*f, ey*f, ez*f))
                if theta - theta_c > p['sat_thr_tilt']:
                    sat_tilt = True

        # ===== 6. Диагностика магнитометра
        s_m = 0.0
        if has_mag:
            mx, my, mz = mag
            nm = sqrt(mx*mx + my*my + mz*mz)
            if nm > 1e-6:
                mh = (mx/nm, my/nm, mz/nm)
                if self.B_ref <= 0.0:
                    self.B_ref = nm
                s_norm_m = self._gauss(nm/self.B_ref - 1.0, p['mag_norm_tol'])
                # дифференциальный тест
                lam = exp(-dt / p['diff_leak_tau'])
                if self._pm is not None:
                    dx_ = mh[0] - self._pm[0] + (wy*mh[2] - wz*mh[1])*dt
                    dy_ = mh[1] - self._pm[1] + (wz*mh[0] - wx*mh[2])*dt
                    dz_ = mh[2] - self._pm[2] + (wx*mh[1] - wy*mh[0])*dt
                    self.cmx = self.cmx*lam + dx_
                    self.cmy = self.cmy*lam + dy_
                    self.cmz = self.cmz*lam + dz_
                cnm = sqrt(self.cmx**2 + self.cmy**2 + self.cmz**2) \
                    / p['diff_leak_tau']
                s_diff_m = self._gauss(cnm, p['mag_diff_tol'] * tol_scale)
                # тест наклонения: â·m̂ (инвариант вращения) против эталона
                if na > 1e-6 and self.gate_a:
                    dip_ref_vec = ah
                else:
                    dip_ref_vec = gravity_body(q)
                dip = (dip_ref_vec[0]*mh[0] + dip_ref_vec[1]*mh[1]
                       + dip_ref_vec[2]*mh[2])
                if not self._dip_init:
                    self.S_dip = dip
                    self._dip_init = True
                s_dip = self._gauss(dip - self.S_dip, p['dip_tol'])
                s_m = s_norm_m * s_diff_m * s_dip
                self._pm = mh
                # онлайн-обучение эталонов поля в доверенные периоды
                if s_m > 0.7 and self.gate_a:
                    k_f = dt / (p['field_tau'] + dt)
                    self.B_ref += (nm - self.B_ref) * k_f
                    self.S_dip += (dip - self.S_dip) * k_f
        else:
            self._pm = None
            self.cmx = self.cmy = self.cmz = 0.0
        self.Em = self._evidence(self.Em, s_m, dt, up_tau=p['ev_up_tau_mag'])
        self.gate_m = self._gate(self.gate_m, self.Em)

        # ===== 7. Коррекция yaw: снап вокруг земной вертикали в пределах бюджета
        sat_yaw = False
        if has_mag and self.gate_m and self.gate_a:
            me = quat_rotate(q, mh)
            h = sqrt(me[0]*me[0] + me[1]*me[1])
            if h > p['mag_hmin']:
                psi = atan2(me[0], me[1])     # ошибка yaw = -psi
                apsi = psi if psi >= 0.0 else -psi
                if apsi > p['deadband']:
                    allow = self.Ty
                    if self.Em < p['burst_evidence']:
                        routine = p['routine_mult'] * rho * dt
                        allow = allow if allow < routine else routine
                    psi_c = apsi if apsi < allow else allow
                    self.Ty -= psi_c
                    sgn = 1.0 if psi >= 0.0 else -1.0
                    q = quat_mult(quat_from_rotvec(0.0, 0.0, sgn*psi_c), q)
                    if apsi - psi_c > p['sat_thr_yaw']:
                        sat_yaw = True

        # ===== 8. Насыщение бюджета => реальный дрейф больше ожидаемого:
        # форсируем неопределённость (петля самонастройки без ковариаций)
        if sat_tilt or sat_yaw:
            self.ub = min(self.ub + p['ub_boost']*dt*10.0, p['ub_max'])

        self.q = quat_normalize(q)
        self.flag_dyn = dyn_flag
        self.flag_mag = (has_mag and not self.gate_m)
        return self.q

    # ------------------------------------------------------------ diagnostics
    def diagnostics(self):
        return dict(Ea=self.Ea, Em=self.Em,
                    bias=(self.bx, self.by, self.bz),
                    Tt=self.Tt, Ty=self.Ty, ub=self.ub,
                    static=self.is_static,
                    gate_a=self.gate_a, gate_m=self.gate_m)
