import numpy as np

def dist (S, T):
    """ The distance between S and T is calculated by sqrt(Trace((S-T)*(S-T))).
    input : nxn symmetric matrices S, T
    output: sqrt(Trace((S-T)*(S-T))) """
    diff = S - T
    return np.sqrt(np.trace(diff.dot (diff)))

def gauss_algorithm (S_input):
    """ Gauss reduction algorithm for 2D lattices.
    input : 2x2 symmetric positive-definite S = S_input
    output: 2x2 basis transform matrix g and g*S*(g^T)
             such that g S g^T = (sij) satisfies 0<=-2*s12<=s11<=s22 """
    ndim = 2
    assert S_input.shape == (ndim,ndim)
    
    S = S_input.copy()
    g = np.identity (ndim, dtype=int) # Set g = I
    arr = np.array ([[0,1], [1,0]])
    while True:
        assert S[1,1] > 0, "The argument S_input is not positive definite."
        # Let m be the integer closest to - s12 / s22
        m = int (np.round (-S[0,1] / S[1,1]))

        # S = ((0,1),(1,m))*S*((0,1),(1,m))
        arr[1][1] = m
        # g = arr*g, arr*S*(arr^T) 
        g = arr.dot (g)
        S = arr.dot (S).dot (arr)
        if S[0,0] <= S[1,1]:
            break
    
    if S[0,1] > 0:
        g[0] *= -1 
        S[0,1] *= -1 
        S[1,0] *= -1
    
    assert 0 < -S[0,1]*2 <= S[0,0] <= S[1,1], "Reduced S:\n" + str(S)
    assert dist(S, g.dot (S_input).dot (g.T)) <= dist(S_input, np.zeros((ndim,ndim)))*1.0e-10
    return g, S


def error_stable_bravais_lattice_detrmination (S_obs):
    """ Error stable bravais lattice determination
    input : 2x2 symmetric positive-definite S = S_obj 
    output: Ans_rP , Ans_hP , Ans_rC, Ans_sP are
            lists of pairs of an integral matrix g
            and the projection of g Sobs tg to a metric tensor 
            with the exact symmetry of the Bravais classes:
            - primitive rectangular (rP), 
            - hexagonal (hP),
            - centered rectangular  (rC),
            - square (sP).                """
    
    assert S_obs.shape == (2,2)

    # (1) Gauss reduction: compute the S := g0*Sobs*(g0^T).
    g0, S = gauss_algorithm (S_obs)

    # (2) Centring: 
    # The following is a list for the centered case.
    hf = np.array ([[1,1],[1,-1]]) # Matrix for centering.
    arr_1 = np.array ([[1,0],[-1,-1]])
    arr_2 = np.array ([[0,1],[-1,-1]])
    # CrC :={ hF, hF*((1 0)(-1 -1)), hF*((0 1)(-1 -1)) }
    CrC = [hf, hf.dot(arr_1), hf.dot(arr_2)]
    
    # (primitive rectangular, rP)
    # 2-3 2-4 Append ( g0, ((s11 0)(0 s22)) ) to AnsrP
    ans_rP = [[ g0, np.diag (np.diag (S)) ]]

    # (hexagonal, Hp)
    # compute s = s11 + s22 + 2 * s12
    # Append (g0, (s,-s/2)(-s/2,s))
    ans_hP = []
    s = S[0,0] + S[1,1] + S[0,1]*2
    ans_hP.append ((g0, np.array ([[s, -s/2],[-s/2, s]])))

    # (centered rectangular, rC)
    ans_rC = []
    for g in CrC:
        # compute Snew = g*S*(g^T)
        Snew = g.dot (S).dot(g.T)
        s11, s22 = np.diag (Snew)

        # put Snew = ((s11,0)(0,s22))
        # Replace s11, s22 in S and the rows of g if s11 > s22
        if s11 > s22:
            s11, s22 = s22, s11
            g = g[::-1, :]
        Snew = np.diag ([s11, s22])
        ans_rC.append ((g.dot(g0), Snew))

    # (3) projection to square (sP)
    ans_sP = []
    for g, S in ans_rP: # Note that ans_rp is length 1.
        # s := (s11 + s22)/2, where si j is the (i; j)-entry of S.
        s = np.sum (np.diag (S)) / 2
        # Append (g, ((s,0)(0,s)))
        ans_sP.append ((g, np.diag ([s, s])))

    return [ans_rP, ans_hP, ans_rC, ans_sP]



if __name__ == '__main__':
    ndim = 2
    # Make a positive-definite matrix S from a basis matrix B.
    B = [np.random.random() for _ in range (ndim*ndim)]
    B = np.array (B).reshape (ndim,ndim)
    S_input = B.dot(B.T)
    print (S_input)

    # error stable bravais lattice determination algorithm
    ans = error_stable_bravais_lattice_detrmination (S_input)
    title = [ "Primitive rectangular", "Hexagonal", "Centered rectangular", "Square" ]
    for i in range(len(ans)):
        print ("* ", title[i], ":")
        for cand in ans[i]:
            g = cand[0]
            Snew = g.dot (S_input).dot(g.T)
            print ("*  g, g*Sinp*(g^T), S_proj, distance = ", dist(Snew, cand[1]))
            for j in range(ndim):
                print ("  ", cand[0][j], Snew[j], cand[1][j])
        print ("")
