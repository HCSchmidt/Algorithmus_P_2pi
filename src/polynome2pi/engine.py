import cmath

class PolynomeEngine:
    """
    Owns all precomputed constants and all scratch state used by energie().
    This eliminates globals while preserving performance (scratch reused).
    """

    def __init__(self, pow_off: int = 40):
        self.pi = cmath.pi
        self.two_pi = 2 * self.pi
        self.pow_off = pow_off

        # Precompute powers (2π)^k into a list for fast indexing
        self.two_pi_pow = [0.0] * (2 * pow_off + 1)
        for k in range(-pow_off, pow_off + 1):
            self.two_pi_pow[k + pow_off] = self.two_pi ** k

        # C-dependent constants
        pi = self.pi
        self.E_C_POS = (
            -pi
            + 2 * pi ** (-1)
            - pi ** (-3)
            + 2 * pi ** (-5)
            - pi ** (-7)
            + pi ** (-9)
            - pi ** (-12)
            - 2 * pi ** (-14)
        )
        self.E_C_NEG = 2 * pi - pi ** (-1) + self.E_C_POS
        self.E_C_ZERO = pi ** (-12) + 2 * pi ** (-14)

        # Scratch state (reused; no per-call allocations)
        self.E = [0.0] * 10
        self.g = [[0.0] * 10 for _ in range(10)]  # no aliasing

        # Hot precomputations for energie()
        off = self.pow_off
        p = self.two_pi_pow

        # (2π)^l for l in {4,3,2}
        self.POW_L_4 = p[4 + off]
        self.POW_L_3 = p[3 + off]
        self.POW_L_2 = p[2 + off]

        # (2π)^n for n in {1,0,-1}
        self.POW_N_1 = p[1 + off]
        self.POW_N_0 = p[0 + off]
        self.POW_N_M1 = p[-1 + off]

        # 2*(2π)^(-8) constant used in E6
        self.POW_NEG8_2 = 2.0 * p[-8 + off]

        # Precompute 9 combos for -l-n-1 and -l-n
        LS = (4, 3, 2)
        NS = (1, 0, -1)

        # include factor 2 for E3/E4/E5 terms
        self.POW_LN_M1 = [[0.0] * 3 for _ in range(3)]      # 2*(2π)^(-l-n-1)
        self.POW_LN_0 = [[0.0] * 3 for _ in range(3)]       # 2*(2π)^(-l-n)

        # exclude factor 2 for E7
        self.POW_LN_M1_ONLY = [[0.0] * 3 for _ in range(3)] # (2π)^(-l-n-1)
        self.POW_LN_0_ONLY = [[0.0] * 3 for _ in range(3)]  # (2π)^(-l-n)

        for li, l in enumerate(LS):
            for ni, n in enumerate(NS):
                self.POW_LN_M1[li][ni] = 2.0 * p[-l - n - 1 + off]
                self.POW_LN_0[li][ni] = 2.0 * p[-l - n + off]
                self.POW_LN_M1_ONLY[li][ni] = p[-l - n - 1 + off]
                self.POW_LN_0_ONLY[li][ni] = p[-l - n + off]

    def energie(self, i4, i3, i2, i1, i0, i_1, C):
        """
        Behavior-identical to the original Energie() nested-loop logic,
        but faster due to:
          - indexed pow tables (no dict)
          - local accumulators, fewer list writes
          - precomputed constants for the 9 (l,n) combos
        """
        g2 = self.g[2]
        g1 = self.g[1]

        # set inputs
        g2[4] = i4
        g2[3] = i3
        g2[2] = i2
        g1[1] = i1
        g1[0] = i0
        g1[-1] = i_1

        # C-dependent base term
        if C > 0:
            E0 = C * self.E_C_POS
        elif C < 0:
            E0 = -C * self.E_C_NEG
        else:
            E0 = self.E_C_ZERO

        # gluons and fermions
        E2 = g2[4] * self.POW_L_4 + g2[3] * self.POW_L_3 + g2[2] * self.POW_L_2
        E1 = -(g1[1] * self.POW_N_1 + g1[0] * self.POW_N_0 + g1[-1] * self.POW_N_M1)

        E3 = 0.0
        E4 = 0.0
        E5 = 0.0
        E6 = 0.0
        E7 = 0.0

        # exact original loop/break semantics
        for li, l in enumerate((4, 3, 2)):
            for ni, n in enumerate((1, 0, -1)):
                gl = g2[l]
                gn = g1[n]

                if gl != 0 and gn != 0:
                    ln = l + n

                    if ln < 4:
                        if gl > 0:
                            E3 += gl * gn * self.POW_LN_M1[li][ni]
                        else:
                            E4 += gl * gn * self.POW_LN_0[li][ni]

                    if ln > 3:
                        E5 -= gl * gn * self.POW_LN_M1[li][ni]

                    prod = gl * gn
                    E6 += (prod if prod >= 0 else -prod) * self.POW_NEG8_2

                    g2[l] = 0
                    g1[n] = 0
                    break

                if gl == 0 and gn == 0:
                    E7 -= self.POW_LN_M1_ONLY[li][ni]
                    E7 -= self.POW_LN_0_ONLY[li][ni]
                    break

        total = E0 + E1 + E2 + E3 + E4 + E5 + E6 + E7

        # preserve side-effects (main loop reads E[0])
        E = self.E
        E[0] = total
        E[1] = E1
        E[2] = E2
        E[3] = E3
        E[4] = E4
        E[5] = E5
        E[6] = E6
        E[7] = E7

        return total