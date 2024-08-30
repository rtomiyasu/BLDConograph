#ifndef _ERROR_STABLE_BRAVAIS_2D_
#define _ERROR_STABLE_BRAVAIS_2D_

# include <iostream>
# include "../utility_data_structure/nrutil_nr.hh"

using namespace std;

template <class T>
void gauss_algorithm (const NRMat<T>& S_input, NRMat<Int4>& g, NRMat<T>& S)
/*Gauss reduction algorithm for 2D lattices.
    input : 2x2 symmetric positive-definite S = S_input
    output: 2x2 basis transform matrix g and g*S*(g^T)
             such that g S g^T = (sij) satisfies 0<=-2*s12<=s11<=s22 */
{
    const Int4 ndim = 2;
    assert (S_input.nrows() == ndim && S_input.ncols() == ndim);

    S = S_input;
    g = identity_matrix<Int4> (ndim);
    NRMat<Int4> arr = put_matrix (vector<vector<Int4>> {{0,1},{1,0}});

    while (true)
    {
        assert (S[1][1] > 0); //"The argument S_input is not positive definite."
        // Let m be the integer closest to - s12 / s22
        Int4 m = round (- S[0][1] / S[1][1]);

        // Let m be the integer closest to - s12 / s22
        arr[1][1] = m;
        //g = arr*g, arr*S*(arr^T) 
        g = mprod (arr, g);
        S = mprod (mprod (arr, S), arr);

        if (S[0][0] <= S[1][1]) break;
    };

    if (S[0][1] > 0)
    {
        g[0][0] *= -1; g[0][1] *= -1;
        S[0][1] *= -1; S[1][0] *= -1;
    };

    assert (0 <= -S[0][1] * 2 && -S[0][1] <= S[0][0] && S[0][0] <= S[1][1]);
    assert (dist (S, mprod (mprod (g, S_input), transpose (g))) <= dist (S_input, NRMat<T>(ndim, ndim, 0))*1.0e-10);
}

#endif
