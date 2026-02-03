import numpy as np

def init_result_arrays():
    D_i_N = [0.0] * 300
    D_i_c_ = [""] * 35
    Cnt = [0] * 520

    Emax = np.zeros((35, 520), dtype=float)
    Emin = np.zeros((35, 520), dtype=float)
    i_Emax = np.zeros((35, 520), dtype=int)
    i_Emin = np.zeros((35, 520), dtype=int)
    Dmax = np.zeros((35, 520, 7), dtype=int)
    Dmin = np.zeros((35, 520, 7), dtype=int)

    return D_i_N, D_i_c_, Cnt, Emax, Emin, i_Emax, i_Emin, Dmax, Dmin


def init_plot_buffers():
    xs_by_j = [[] for _ in range(27)]
    ys_by_j = [[] for _ in range(27)]
    grey_segments = []
    return xs_by_j, ys_by_j, grey_segments