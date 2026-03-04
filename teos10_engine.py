# teos10_general_engine.py
import numpy as np
import re
import math
import gsw
from pathlib import Path

class TEOS10_General_Engine:
    """
    Python 严格版 TEOS10 高阶通用引擎
    完全等价于你的 MATLAB Final Strict Version
    """

    def __init__(self):
        self.scales = {
            "s": 0.0248826675584615,
            "t": 0.025,
            "p": 1e-4
        }
        self._load_gsw_coefficients()

    # ==============================
    # 公共接口
    # ==============================
    def calculate_all_orders(self, SA, CT, p, max_order=8):
        dT = self._calc_leibniz(SA, CT, p, max_order, "CT")
        dS = self._calc_leibniz(SA, CT, p, max_order, "SA")

        out = {}
        for n in range(1, max_order + 1):
            out[f"d{n}_T"] = dT[:, n]
            out[f"d{n}_S"] = dS[:, n]
        return out

    # ==============================
    # Leibniz 递推
    # ==============================
    def _calc_leibniz(self, SA, CT, p, N, var):
        SA = SA.reshape(-1)
        CT = CT.reshape(-1)
        p  = p.reshape(-1)

        v_derivs = self._eval_v_derivs(SA, CT, p, N, var)

        num = len(SA)
        rho_derivs = np.zeros((num, N + 1))
        v0 = v_derivs[:, 0]
        rho0 = 1.0 / v0
        rho_derivs[:, 0] = rho0

        for n in range(1, N + 1):
            acc = np.zeros(num)
            for k in range(n):
                acc += math.comb(n, k) * rho_derivs[:, k] * v_derivs[:, n - k]
            rho_derivs[:, n] = -rho0 * acc

        return rho_derivs

    # ==============================
    # 幂函数求导
    # ==============================
    def _eval_v_derivs(self, SA, CT, p, N, target):
        num = len(SA)
        v = np.zeros((num, N + 1))

        s_fac, t_fac, p_fac = self.scales.values()
        Xs = SA + 24
        Xt = CT
        Xp = p

        for order in range(N + 1):
            acc = np.zeros(num)
            for (i, j, k, val) in self.coeffs:
                if target == "CT":
                    if j < order:
                        continue
                    coeff = math.prod(range(j, j - order, -1)) if order > 0 else 1
                    term = (
                        val *
                        (s_fac ** (i / 2)) * (Xs ** (i / 2)) *
                        (p_fac ** k) * (Xp ** k) *
                        (t_fac ** j) * coeff * (Xt ** (j - order))
                    )

                else:  # SA
                    exp = i / 2
                    if exp < order:
                        continue
                    coeff = math.prod([exp - q for q in range(order)]) if order > 0 else 1
                    term = (
                        val *
                        (t_fac * Xt) ** j *
                        (p_fac * Xp) ** k *
                        (s_fac ** (i / 2)) *
                        coeff * (Xs ** (exp - order))
                    )

                acc += term

            v[:, order] = acc

        return v

    # ==============================
    # 系数读取（关键）
    # ==============================
    def _load_gsw_coefficients(self):
        import inspect
        import gsw

        path = Path(inspect.getfile(gsw)) / "specvol.py"
        text = path.read_text()

        rgx = re.compile(
            r"^\s*v(\d)(\d)(\d)\s*=\s*([-+eE\d\.]+);",
            re.MULTILINE
        )

        self.coeffs = np.array([
            (int(a), int(b), int(c), float(d))
            for a, b, c, d in rgx.findall(text)
        ])

        if self.coeffs.shape[0] != 48:
            print(f"⚠️ 系数数量异常: {self.coeffs.shape[0]}（应为48）")
