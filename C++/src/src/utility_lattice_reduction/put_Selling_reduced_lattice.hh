/*
 * The MIT License

   Conograph (powder auto-indexing program)

Copyright (c) <2012> <Ryoko Oishi-Tomiyasu, KEK>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 *
 */
#ifndef _PUT_SELLING_REDUCED_LATTICE_HH_
#define _PUT_SELLING_REDUCED_LATTICE_HH_

#include <iostream>
#include "matrix_NbyN.hh"
#include "../utility_data_structure/index_set.hh"
#include "../utility_data_structure/nrutil_nr.hh"
#include "../utility_func/transform_sym_matrix.hh"
#include "../bravais_type/matrix_3by3.hh" // put_perm_mat3

using namespace std;

inline Double dist (const NRMat<Double>& Sp, const NRMat<Double>& Ta)
/*The distance between S and T is calculated by sqrt(Trace((S-T)*(S-T))).
    input : nxn symmetric matrices S, T
    output: sqrt(Trace((S-T)*(S-T)))*/
{
    const NRMat<Double> diff = difference (Sp, Ta); // S - T
    return sqrt (max (0.0, (trace (mprod (diff, diff)))));
}

inline bool equiv_eps (const Double& s, const Double& t, const Double& eps)
/* Check whether s and t are nearly equal, i.e., abs(s-t) <= eps*min(s,t).
    s, t  : values,
    eps   : threshold on relative error,
    output: True or False.*/
{
    const Double max_ST = max (s, t);
    const Double diff_ST = abs (s - t);
    return diff_ST <= eps * max_ST;
}

inline bool check_equiv (const NRMat<Double>& Sp, const NRMat<Double>& Ta, const Double& eps)
{
    const Int4 irow = Sp.nrows();
    const Int4 icol = Sp.ncols();
    for (Int4 i = 0; i < irow; i++)
    {
        if ( !equiv_eps (Sp[i][i], Ta[i][i], eps)) return false;
        for (Int4 j = i + 1; j < icol; j++)
        {
            if (!equiv_eps (
            Sp[i][i] + Sp[j][j] + 2. * Sp[i][j],
            Ta[i][i] + Ta[j][j] + 2. * Ta[i][j], eps)) return false;
        }
    }
    return true;
}

inline NRMat<Int4> I_ext ()
// [[1,0,0],[0,1,0],[0,0,1],[-1,-1,-1]]
{
    Int4 irow = 4; Int4 icol = 3;
    NRMat<Int4> ans (irow, icol, 0);
    ans[0][0] = ans[1][1] = ans[2][2] = 1;
    ans[3][0] = ans[3][1] = ans[3][2] =-1;
    return ans;
}

template <class T>
NRMat<T> S_3x3_to_4x4 (const NRMat<T>& S)
{
    NRMat<Int4> In = I_ext ();
    return mprod (mprod (In, S), transpose (In));
}

template <class T>
bool positive_index (const NRMat<T>& S, vector<Int4>& idx)
{
    bool flg = false;
    Int4 i; Int4 j;
    for (i = 0; i < S.nrows(); i++)
    {
        for (j = i + 1; j < S.ncols(); j++)
        {
            if (S[i][j] > 0)
            {
                flg = true;
                break;
            }
        }
    if (flg) break;
    }

    idx = {i, j};
    return flg;
}


inline bool Selling_reduction_3D (const NRMat<Double>& S_input, NRMat<Int4>& g, NRMat<Double>& S)
/*Selling reduction algorithm for 3D lattices from Conway's book.
    input : 4x4 matrix S = S_input equal to I_ext*S*I_ext^T of a symmetric positive-definite S. 
    output: 4x4 basis transform matrix g and g*S*(g^T)
             such that g S g^T is Selling reduced (i.e., all non-diagonal entries are <= 0).*/
{
    Int4 ISIZE = 4;
    assert ((S_input.nrows() == ISIZE) & 
            (S_input.ncols() == ISIZE));
    S = S_input;
    g = identity_matrix<Int4> (ISIZE); // set g = I
    Int4 itnum = 0;
    while (true)
    {
        if ((S[0][0] <= 0 || S[1][1] <= 0 || S[2][2] <= 0) || (itnum > 10000))
        {
            return false;
        };
        itnum += 1;

        // Search for i, j such that S[i,j] > 0
        vector<Int4> idx; bool flg; Int4 i, j;
        flg = positive_index (S, idx);
        i = idx[0]; j = idx[1];
        //# If non-diagonal entries of S are <= 0,
        // Selling reduction is successful.
        if (!flg)
        {
            NRMat<Double> zeros (ISIZE, ISIZE, 0.0); 
            Double d1 = dist (S, mprod (mprod (g, S_input), transpose (g)));
            Double d2 = dist (S_input, zeros);
            assert (d1 <= d2);
            return true;
        }
        else
        {
            NRMat <Int4> g2 = put_reduction_matrix_4D (i, j);
            g = mprod (g2, g);
            S = mprod (mprod (g2, S), transpose (g2));
        }
    };
    return false;
}

template <class T>
NRMat<Int4> putMatrixToMoveSmallerDiagonalLeftUpper (const NRMat<T>& S)
{
    Int4 ndim = S.nrows(); 

    // take diagonal entries of S
    vector<T> diag (ndim);
    for (Int4 i = 0; i < ndim; i++) diag[i] = S[i][i];

    // Sort in ascending order
    vector<Int4> index = argsort (diag);
    NRMat <Int4> g (ndim, ndim, 0);
    for (Int4 k = 0; k < ndim; k++)
    {
        g[k][index[k]] = 1;
    }
    return g;
}

template<class T>
void moveSmallerDiagonalLeftUpper (NRMat<Int4>& g, NRMat<T>& S)
/*Move smaller diagonal to left upper direction
    input: nxn matrix S
    output: permutation matrix g and g*S*(g^T)
    such that the entries of g S g^T are in ascending order. */
{
    g = putMatrixToMoveSmallerDiagonalLeftUpper (S);
    S = mprod (mprod (g, S), transpose (g));
}

template <class T>
void Delaunay_reduction (
        const NRMat<T>& S_input, NRMat<Int4>& g, NRMat<T>& S)
/*Delaunay reduction algorithm for 3D lattices.
    input : 3x3 symmetric positive-definite matrix S = S_input. 
    output: 3x3 basis transform matrix g and g*S*(g^T)
             such that S_red = g S g^T = (sij) is Delaunay reduced (i.e., s11<=s22<=s33<=s44
             and all non-diagonal entries of I_ext*S_red*I_ext^T are <= 0). */
{
    NRMat<T> S4_ = S_3x3_to_4x4 (S_input);
    bool flg; NRMat<Int4> g_4x4; NRMat<Double> S4;
    flg = Selling_reduction_3D (S4_, g_4x4, S4);
    assert (flg); // Selling reduction should not be failed.
    NRMat<Int4> g2_4x4;
    moveSmallerDiagonalLeftUpper (g2_4x4, S4);

    // Remove the last column of g_44*I_ext and
    // the last row & column of S4.
    g = mprod (mprod (g2_4x4, g_4x4), I_ext());
    g = mat_delete (g, 3, 0);
    S = mat_delete (mat_delete (S4, 3, 0), 3, 1);
    NRMat<T> zeros (3, 3, 0.0);

    Double d1 = dist (S, mprod (mprod (g, S_input), transpose (g)));
    Double d2 = dist (S_input, zeros);
    assert  (d1 <= d2);
    //return forward_as_tuple (g, S);
}

template<class T>
void Delaunay_reduction_of_inverse (
                const NRMat<T>& S, NRMat<Int4>& g, NRMat<Double>& S_del_inv)
/*Selling reduction of S^{-1}.
    input : 3x3 symmetric positive-definite S. 
    output: 3x3 basis transform matrix g 
    such that (g S g^T)^{-1} is Selling reduced
    (i.e., all non-diagonal entries are <= 0). */
{
    NRMat<T> Sinv = inverse_mat_3x3 (S);
    NRMat<Int4> h; NRMat<Double> Sinv_red;
    Delaunay_reduction (Sinv, h, Sinv_red);
    g = inverse_mat_3x3 (transpose (h));
    S_del_inv = mprod (mprod (g, S), transpose (g));
    //return g, S_del_inv
}


#endif /*SUPER_BASIS3_HH_*/
