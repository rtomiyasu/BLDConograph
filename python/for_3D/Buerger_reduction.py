import numpy as np

def dist (S, T):
    """ The distance between S and T is calculated by sqrt(Trace((S-T)*(S-T))).
    input : nxn symmetric matrices S, T
    output: sqrt(Trace((S-T)*(S-T))) """
    diff = S - T
    return np.sqrt(max(0.,np.trace(diff.dot (diff))))

def equiv_eps (s, t, eps):
    """ Check whether s and t are nearly equal, i.e., abs(s-t) <= eps*min(s,t).
        s, t  : values,
        eps   : threshold on relative error,
        output: True or False. """
    min_ST = np.minimum (s, t)
    diff_ST = np.abs (s - t)
    return diff_ST <= eps * min_ST

def check_equiv (S, T, eps):
    """ Check whether S and T are nearly equal, i.e., for every 1 <= i < j <= n,
        (i) |sii - tii| <= eps*min (sii, tii)
        (ii)|(sii + sjj + 2sij) - (tii + tjj + 2tij)| <= eps*min (sii + sjj + 2sij, tii + tjj + 2tij)
        S, T  : nxn symmetric matrices,
        eps   : threshold on relative error,
        output: True or False. """
    assert S.shape == T.shape
    nrow, ncol = S.shape
    for i in range(nrow):
        if not equiv_eps(S[i,i], T[i,i], eps):
            return False
        for j in range(i+1, nrow):
            if not equiv_eps(S[i,i] + S[j,j] + 2. * S[i,j], T[i,i] + T[j,j] + 2. * T[i,j], eps):
                return False
    return True


def Buerger_reduction (S_input):
    """ Buerger reduction (= 3D case of Minkowski reduction) algorithm for 3D lattices.
    input : 3x3 symmetric positive-definite S = S_input
    output: 3x3 basis transform matrix g and g*S*(g^T)
             such that g S g^T = (sij) satisfies 
             - s11<=s22<=s33,
             - 2*|s12|,2*|s13| <= s11, 2*|s23| <= s22,
             - s12, s13, s23 > 0 or s12, s13, s23 <= 0,
             - -2*(s12 + s13 + s23) <= s11+s22. """
    ndim = 3
    arr_list = [ np.array ([[0,1,0], [1,0,0], [0,0,1]]),
                 np.array ([[0,0,1], [0,1,0], [1,0,0]]),
                 np.array ([[1,0,0], [0,0,1], [0,1,0]]),
                 np.identity (ndim, dtype=int) ] # these matrices are all symmetric.
    arr_list[3][2,0] = 1
    arr_list[3][2,1] = 1
    # [[1,0,0],[0,1,0],[1,1,1]]

    assert S_input.shape == (ndim,ndim)
    S = S_input.copy()
    g = np.identity (ndim, dtype=int) # Set g = I
    while True:
        k=0
        # 0≦i＜j≦3
        for i in range(ndim): # [0,1,2]
            for j in range(i+1,ndim): # (i,j) = (0,1),(0,2),(1,2)
                arr = arr_list[k].copy(); k += 1
                if S[i,i] < S[j,j]:
                    g = arr.dot (g)
                    S = arr.dot (S).dot (arr)
                # Let m be the integer closest to - sij / sjj
                assert S[j,j] > 0, "The argument S_input is not positive definite:\n "+str(S)
                m = int (np.round (- S[i,j] / S[j,j]))
                # E.g., if (i,j)=(0,1), arr = ((0,1,0),(1,m,0),(0,0,1))*S*((0,1,0),(1,m,0),(0,0,1))
                arr[j,j] = m
                # g = arr*g, arr*S*(arr^T) 
                g = arr.dot (g)
                S = arr.dot (S).dot (arr)
        if ( S[0,0] <= S[1,1] <= S[2,2] and 2*abs(S[0,1]) <= S[0,0] and
            2*abs(S[0,2]) <= S[0,0] ):
            assert 2*abs(S[1,2]) <= S[1,1]
            # Make s12, s23 non-positive. 
            diag = [1,1,1]
            if S[0,1] > 0:
                diag[0] = -1
            if S[1,2] > 0:
                diag[2] = -1
            arr = np.diag(diag)
            g = arr.dot (g)
            S = arr.dot (S).dot (arr)
            if 2*(S[0,1] + S[0,2] + S[1,2]) + S[0,0] + S[1,1] < 0: 
                arr = arr_list[3]
                g = arr.dot (g)
                S = arr.dot (S).dot (arr.T)
            else:
                break
                
    # At this point, s12, s23 are non-positive. 
    # If s13 is positive, make s12, s23 positive, unless s12 or s23 = 0.
    diag = [1,1,1]
    if S[0,2] > 0:
        if S[0,1] == 0:
            diag[0] = -1
        elif S[1,2] == 0:
            diag[2] = -1
        else:
            diag[0] = -1
            diag[2] = -1
    arr = np.diag(diag)
    g = arr.dot (g)
    S = arr.dot (S).dot (arr)
    assert ( S[0,0] <= S[1,1] <= S[2,2] and 
            2*abs(S[0,1]) <= S[0,0] and
            2*abs(S[0,2]) <= S[0,0] and
            2*abs(S[1,2]) <= S[1,1] and
            ( (S[0,1] > 0 and S[0,2] > 0 and S[1,2] > 0) or
              (S[0,1]<= 0 and S[0,2]<= 0 and S[1,2]<= 0) ) ), "Reduced S:\n" + str(S)
    assert 2*(S[0,1] + S[0,2] + S[1,2]) + S[0,0] + S[1,1] >= 0, "Reduced S:\n" + str(S)
    assert dist(S, g.dot (S_input).dot (g.T)) <= dist(S_input, np.zeros((ndim,ndim)))*1.0e-10, "Reduced S:\n" + str(S)
    return g, S
    

if __name__ == '__main__':
    ndim = 3
    # Make a positive-definite matrix S from a basis matrix B.
    B = [np.random.random() for _ in range (ndim*ndim)]
    B = np.array (B).reshape(ndim,ndim)
    S_input = B.dot (B.T)
    
    C = [np.random.random() for _ in range (ndim*ndim)]
    C = np.array (C).reshape(ndim,ndim)
    T = C.dot (C.T)

    print (S_input)
    print (T)
    flg = check_equiv (S_input, T, 0.1)
    print (flg)
    """print ("* Input S:")
    print (S_input)

    g, S = Buerger_reduction (S_input)
    print ("* g, Buerger reduced S = g*S*g^T")
    for j in range(ndim):
         print (" ", g[j], S[j])"""
