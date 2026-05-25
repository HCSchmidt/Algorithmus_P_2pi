import cmath
from .cli import ScanSector, PolynomeConfig
def run_scan(
    *,
    config: PolynomeConfig,
    sector: ScanSector,
    obj_E,
    obj_min,
    obj_max,
    Cnt,
    Emax,
    Emin,
    i_Emax,
    i_Emin,
    Dmax,
    Dmin,
    PEmax,
    PEmin,
    xs_by_j,
    ys_by_j,
    grey_segments,
    N_segments,
    data,
):
    
    """
    Executes the nested scan loops. Returns (i_T, i_T1, mmax).
    Mutates the provided arrays/buffers.
    """
    engine = PolynomeEngine()
    energie = engine.energie
    E_local = engine.E
    
    is_heavy = (sector is ScanSector.heavy)
    is_nucleon = (sector is ScanSector.nucleon)
    is_E112P = (sector is ScanSector.E112P)
    is_E333U = (sector is ScanSector.E333U)
    is_E333D = (sector is ScanSector.E333D)
    is_E222 = (sector is ScanSector.E222)

    PE = ""
    i_T = 0
    i_T1 = 0
    mmax = 0
    ct = 0
    pi = cmath.pi ; E1700= (2 * pi)**4 +(2 * pi)**2+ (2 * pi)**1; 

    for i4 in range(-2 * config.J4, 2 * config.J4 + 1):
        i4h = 0.5 * i4
        for i3 in range(-2 * config.J3, 2 * config.J3 + 1):
            i3h = 0.5 * i3
            for i2 in range (-2 * config.J2, 2 * config.J2 + 1):    #(-4,5):  #
                i2h = 0.5 * i2

                for i1 in range (-6, 7):#
                    i1h = 0.5 * i1
                    for i0 in range(-6, 7):   #
                        i0h = 0.5 * i0
                        for i_1 in range(-6, 7):   #
                            i_1h = 0.5 * i_1
                            for C in range(-2, 3):
                                Ch = 0.5 * C 
                                PE = energie(i4h, i3h, i2h, i1h, i0h, i_1h, Ch)
                                E0 = E_local[0]
                                if E0 < 0: continue
                                m = int(512 + 64 * i4 + 8 * i3 + i2)
                                if m > mmax:
                                    mmax = m
                                    ct = 0
                                ct += 1
                                Cnt[mmax] = ct

                                if is_E333U and E0 < 1700: continue
                                if is_E112P and E0 < E1700: continue                                
                                if is_E333D and E0 > 2000: continue 
                                if is_E222 and E0 > 3000: continue 
                                if is_nucleon and (E0 < 1830 or E0 > 1843): continue
                                
                                i_T += 1
                                N_= 10*(abs(i4h)*16 + abs(i3h)*8 + abs(i2h)*4 + abs(i1h)*2 + abs(i0h) + abs(i_1h)/2)
                                #N_= 20*(abs(i4h) + abs(i3h) + abs(i2h) + abs(i1h) + abs(i0h) + abs(i_1h))
                                N_segments.append([(i_T, N_), (i_T + 1, N_)])
                                flag_match = 0

                                for j in range(1, 29):
                                    Ej = obj_E[j]
                                    if (E0 - Ej <= obj_max[j]) and (E0 - Ej >= obj_min[j]):
                                        i_T1 += 1

                                        if Emax[j, m] <= E0:
                                            Emax[j, m] = E0
                                            i_Emax[j, m] = i_T
                                            Dmax[j, m, 0] = i4
                                            Dmax[j, m, 1] = i3
                                            Dmax[j, m, 2] = i2
                                            Dmax[j, m, 3] = i1
                                            Dmax[j, m, 4] = i0
                                            Dmax[j, m, 5] = i_1
                                            Dmax[j, m, 6] = C
                                            PEmax[j][m] = PE

                                        if Emin[j, m] >= E0 or Emin[j, m] == 0:
                                            Emin[j, m] = E0
                                            i_Emin[j, m] = i_T
                                            Dmin[j, m, 0] = i4
                                            Dmin[j, m, 1] = i3
                                            Dmin[j, m, 2] = i2
                                            Dmin[j, m, 3] = i1
                                            Dmin[j, m, 4] = i0
                                            Dmin[j, m, 5] = i_1
                                            Dmin[j, m, 6] = C
                                            PEmin[j][m] = PE

                                        xs_by_j[j].append(i_T)
                                        ys_by_j[j].append(E0)
                                        flag_match = 1
                                        data.append([i_T, E0 , j])
                                        

                                if E0 > 0 and flag_match == 0:
                                    grey_segments.append([(i_T, E0), (i_T + 1, E0)])
                                    data.append([i_T, E0 , 0])
                                    
    return i_T, i_T1, mmax
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
        self.POW_L_4 = p[4 + off]; 
        self.POW_L_3 = p[3 + off]; 
        self.POW_L_2 = p[2 + off]; 
        
        # (2π)^n for n in {1,0,-1}
        self.POW_N_1 = p[1 + off]; 
        self.POW_N_0 = p[0 + off]; 
        self.POW_N_M1 = p[-1 + off]; 

        # 2*(2π)^(-8) constant used in E6
        self.POW_NEG8_2 = 2.0 * p[-8 + off]; 

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
        PE0_POS = " (- π + 2π^-1 - π^-3 + 2π^-5 - π^-7 + π^-9 - π^-12 - 2π^-14) "
        PE0_NEG = " (π + π^-1 - π^-3 + 2π^-5 - π^-7 + π^-9 - π^-12 - 2π^-14) "
        PE0_ZERO = " (+ π^-12 + 2π^-14)"
        
        if C > 0:
            E0 = C * self.E_C_POS;  PE0 = " + " + str(C) + PE0_POS
        elif C < 0:
            E0 = -C * self.E_C_NEG; PE0 = " + " + str(-C) + PE0_NEG
        else:
            E0 = self.E_C_ZERO;  PE0 = " + " + PE0_ZERO
        
        # gluons and fermions
        E1 = g2[4] * self.POW_L_4 + g2[3] * self.POW_L_3 + g2[2] * self.POW_L_2;
        E2 = -(g1[1] * self.POW_N_1 + g1[0] * self.POW_N_0 + g1[-1] * self.POW_N_M1)
        PE3=""; PE4=""; PE5=""; PE6=""; PE7=""
        PE1 = f"{g2[4]} (2π)^4 + {g2[3]}(2π)^3 + {g2[2]}(2π)^2 "
        PE2 = f" -({g1[1]} (2π)^1 + {g1[0]}(2π)^0 + {g1[-1]}(2π)^-1) "
        E3 = 0.0
        E4 = 0.0
        E5 = 0.0
        E6 = 0.0
        E7 = 0.0

#        if i4 == 0 and i3 == 1 and i3 == 0:  print(E0, i4, i3, i2, i1, i0, i_1) 

        # exact original loop/break semantics
        for li, l in enumerate((4, 3, 2)):
            for ni, n in enumerate((1, 0, -1)):
                gl = g2[l]
                gn = g1[n]

                if gl != 0 and gn != 0:
                    ln = l + n

                    if ln < 4:
                        if gl > 0:
                            E3 += gl * gn * self.POW_LN_M1[li][ni];  PE3 += f" + {gl*gn}*2(2π)^({-li-ni-1}) "  # PE3 += " + " + str(gl*gn) + "* 2(2π)^"+str(-li-ni-1)
                        else:
                            E4 += gl * gn * self.POW_LN_0[li][ni];   PE4 += f" + {gl*gn}*2(2π)^({-li-ni}) "

                    if ln > 3:
                        E5 -= gl * gn * self.POW_LN_M1[li][ni];      PE5 += f" - {gl*gn}*2(2π)^({-li-ni-1}) "  # PE5 += " + " + str(-gl*gn) + "* 2(2π)^"+str(-li-ni-1)

                    prod = gl * gn
                    E6 += (prod if prod >= 0 else -prod) * self.POW_NEG8_2; PE6 += f" + {abs(gl*gn)}*2(2π)^(-8) "  #PE6 += " + " + str(abs(prod)) + "* 2(2π)^-8"
 
                    g2[l] = 0
                    g1[n] = 0
                    break

                if gl == 0 and gn == 0:
                    E7 -= self.POW_LN_M1_ONLY[li][ni];   PE7 += f" - (2π)^({-li-ni-1}) "  
                    E7 -= self.POW_LN_0_ONLY[li][ni];    PE7 += f" - (2π)^({-li-ni}) " 
                    break

        total = E0 + E1 + E2 + E3 + E4 + E5 + E6 + E7;   PE = PE1 + PE2 + PE3 + PE4 + PE5 + PE6 + PE7 + PE0

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

        return PE 
