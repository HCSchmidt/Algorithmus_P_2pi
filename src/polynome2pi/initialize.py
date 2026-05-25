import numpy as np

def init_result_arrays():
    D_i_N = [0.0] * 300                         #alt 300
    D_i_c_ = [""] * 35
    Cnt = [0] * 1030                           # alt 520

    Emax = np.zeros((35, 1030), dtype=float)
    Emin = np.zeros((35, 1030), dtype=float)     # alle alt 520
    i_Emax = np.zeros((35, 1030), dtype=int)
    i_Emin = np.zeros((35, 1030), dtype=int)
    Dmax = np.zeros((35, 1030, 7), dtype=int)
    Dmin = np.zeros((35, 1030, 7), dtype=int)
    PEmax = [[""] * 1030] * 35
    PEmin = [[""] * 1030] * 35
    return D_i_N, D_i_c_, Cnt, Emax, Emin, i_Emax, i_Emin, Dmax, Dmin, PEmax, PEmin

def init_plot_buffers():
    xs_by_j = [[] for _ in range(29)]
    ys_by_j = [[] for _ in range(29)]
    grey_segments = []
    N_segments = []
    data = []
    return xs_by_j, ys_by_j, grey_segments, N_segments, data
