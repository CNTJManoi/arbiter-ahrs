"""
vqf_wrapper.py — обёртка официальной реализации VQF (Laidig & Seel 2023)
под общий API проекта. Используется ТОЛЬКО на BROAD (равномерная
дискретизация): официальная реализация VQF предполагает фиксированный Ts,
синтетический сценарий проекта — с неравномерным dt по построению.
"""
import numpy as np
from vqf import VQF


class VQFWrapper:
    name = 'VQF'

    def __init__(self, tauAcc=None, tauMag=None):
        self._tauAcc = tauAcc
        self._tauMag = tauMag
        self._f = None
        self.q = (1.0, 0.0, 0.0, 0.0)

    def update(self, gyr, acc, mag, dt):
        if self._f is None:
            kw = {}
            if self._tauAcc is not None:
                kw['tauAcc'] = self._tauAcc
            if self._tauMag is not None:
                kw['tauMag'] = self._tauMag
            self._f = VQF(gyrTs=float(dt), **kw)
        g = np.asarray(gyr, dtype=float)
        a = np.asarray(acc, dtype=float)
        m = np.asarray(mag, dtype=float)
        if np.all(np.isfinite(m)):
            self._f.update(g, a, m)
            q = self._f.getQuat9D()
        else:
            self._f.update(g, a)
            q = self._f.getQuat9D()
        self.q = (q[0], q[1], q[2], q[3])
        return self.q
