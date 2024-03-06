import numpy as np
from Selling_Delaunay_reducion import Delaunay_reduction

def Delaunay_to_Dirichlet (S_del):
    """ Transform Delaunay reduced S_del to a Dirichlet reduced symmetric matrix
        by the method in Section 4 of Balashov & Ursell (1957).
    input : Delaunay reduced S_del (3x3 symmetric positive-definite matrix)
    output: 3x3 basis transform matrix g and g*S*(g^T)
             such that g S g^T = (sij) is Dirichlet reduced. """
    g = np.identity (3, dtype=int) # Set g = I
    if   S_del[0,0] + S_del[0,1]*2 < 0: # Case of A + 2H < 0
        g[1,0] = 1
    elif S_del[0,0] + S_del[0,2]*2 < 0: # Case of A + 2G < 0
        g[2,0] = 1
    elif S_del[1,1] + S_del[1,2]*2 < 0: # Case of B + 2F < 0
        g[2,1] = 1
    S = g.dot (S_del).dot (g.T)
    assert ( S[0,0] <= S[1,1] <= S[2,2] and 
            2*abs(S[0,1]) <= S[0,0] and
            2*abs(S[0,2]) <= S[0,0] and
            2*abs(S[1,2]) <= S[1,1] and
            2*abs(S[0,1] + S[0,2] + S[1,2]) <= S[0,0] + S[1,1] ), "Reduced S:\n" + str(S)
    return g, S


def Delaunay_to_Buerger (S_del):
    """ Transform Delaunay reduced S_del to a Buerger reduced symmetric matrix  
        using the above Delaunay_to_Dirichlet.
    input : Delaunay reduced S_del (3x3 symmetric positive-definite matrix)
    output: 3x3 basis transform matrix g and g*S*(g^T)
             such that g S g^T = (sij) is Buerger reduced. """
    g, S = Delaunay_to_Dirichlet (S_del)
    # Make s12, s23 non-positive, unless s13 becomes positive. 
    # If this makes s13 positive, make s12, s23 positive, unless s12 or s23 = 0.
    arr = np.identity (3, dtype=int) # Set g = I
    s13 = S[0,2]
    if S[0,1] > 0:
        arr[0] = -1; s13 *= -1
    if S[1,2] > 0:
        arr[2] = -1; s13 *= -1
    if s13 > 0:
        if S[0,1] != 0:
            arr[0] *= -1
        if S[1,2] != 0:
            arr[2] *= -1
    g = arr.dot (g)
    S = arr.dot (S).dot (arr)
    assert ( (S[0,1] > 0 and S[0,2] > 0 and S[1,2] > 0) or
              (S[0,1]>= 0 and S[0,2]>= 0 and S[1,2]>= 0) ), "Reduced S:\n" + str(S)
    return g, S

if __name__ == '__main__':
    ndim = 3
    # Make a positive-definite matrix S from a basis matrix B.
    B = [np.random.random() for _ in range (ndim*ndim)]
    B = np.array (B).reshape(ndim,ndim)
    S_input = B.dot (B.T)
    print ("* Input S:")
    print (S_input)

    g, S = Delaunay_reduction (S_input)
    g2, S = Delaunay_to_Buerger (S)
    g = g2.dot (g)
    print ("* g, Buerger reduced S = g*S*g^T")
    for j in range(ndim):
         print (" ", g[j], S[j])
