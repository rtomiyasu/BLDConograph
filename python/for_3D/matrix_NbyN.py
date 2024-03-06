import numpy as np

def put_complement_set(dim : int, i : int, j : int):
    assert i != j
    assert (0 <= i < dim) and (0 <= j < dim)
    flag_tray = [False for _ in range (dim)]
    flag_tray[i], flag_tray[j] = True, True
    indices = np.arange (dim)
    return indices[(indices != i) & (indices != j)]

def put_reduction_matrix_4D (i : int, j : int):
    indices = put_complement_set(4, i, j)
    k, l = indices[0], indices[1]
    transmat = np.zeros ((4,4), dtype = int)
    transmat[i,i], transmat[j,j], transmat[k,k], transmat[l,l] = -1, 1, 1, 1
    transmat[k,i] = 1
    transmat[l,i] = 1
    return transmat

