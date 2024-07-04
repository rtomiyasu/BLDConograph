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
#ifndef PUT_Buerger_REDUCED_LATTICE_HH_
#define PUT_Buerger_REDUCED_LATTICE_HH_

#include <iostream>
#include "../utility_data_structure/SymMat.hh"
#include "../utility_data_structure/nrutil_nr.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../bravais_type/matrix_3by3.hh" // put_perm_mat3

using namespace std;

template<class T>
string mat2text (const NRMat<T>& mat)
{
    string text = "";
    Int4 irow = mat.nrows();
    Int4 icol = mat.ncols();
    for (Int4 i = 0; i < irow; i++)
    {
        for (Int4 j = 0; j < icol; j++)
        {
            text += to_string (mat[i][j]);
            if (j < icol - 1) text += ", ";
            else if (i < irow - 1) text += "\n";
                };
        }
    
    return text;
}

template<class T>
void print_mat (const NRMat<T>& mat)
{
    string text = mat2text (mat);
    cout << text << endl; 
    cout << "nrow : " << mat.nrows() << " ncol : " << mat.ncols() << endl;
}

template<class T>
void buerger_reduction (
            const NRMat<T>& S_input, NRMat<Int4>& g, NRMat<Double>& S)
{
    Int4 ndim = 3;
    //NRMat<Int4> identity_mat = identity_matrix<Int4> (ndim);
    vector<NRMat<Int4>> array_list = {
        put_perm_mat3 (1,0,2), put_perm_mat3 (2,1,0),
        put_perm_mat3 (0,2,1), put_perm_mat3 (0,1,2)
                                                    };
    array_list[3][2][0] = 1;
    array_list[3][2][1] = 1;

    S = S_input;
    g = identity_matrix<Int4> (ndim);

    static const double EPS = 1.0+1.0e-14;
    while (true)
    //for (int l = 0; l < 3; l++)
    {
        Int4 k = 0;
        for (Int4 i = 0; i < ndim; i++) // 0,1,2
        {
            for (Int4 j = i + 1; j < ndim; j++) // 1,2
            // i = 2, skip
            {
                NRMat<Int4> arr = array_list [k];
                k += 1;
                if (S[i][i] < S[j][j])
                {
                    g = mprod (arr, g);
                    S = mprod (mprod (arr, S), arr);
                }
                assert (S[j][j] > 0);
                Int4 m = round (- S[i][j] / S[j][j]);
                arr[j][j] = m;
                // g = arr @ g, S = arr @ S @ arr (<-- arr.T = arr)
                g = mprod (arr, g);
                S = mprod (mprod (arr, S), arr);
            }
        };
        
        if ((S[0][0] <= S[1][1] && S[1][1] <= S[2][2]) & (2 * abs (S[0][1]) <= S[0][0]*EPS) & (2 * abs (S[0][2]) <= S[0][0]*EPS))
            {
                assert (2 * abs (S[1][2]) <= S[1][1]);
                vector<Int4> diag = {1, 1, 1};
                if (S[0][1] > 0) diag[0] = -1;
                if (S[1][2] > 0) diag[2] = -1;
                NRMat<Int4> arr = diagonal (diag);
                g = mprod (arr, g);
                S = mprod (mprod (arr, S), arr);
        
                if( -2 * (S[0][1] + S[0][2] + S[1][2]) <= (S[0][0] + S[1][1])*EPS )
                {
			break;
                }
                else
                { 
                    arr = array_list[3];
                    g = mprod (arr, g);
                    S = mprod (mprod (arr, S), transpose (arr));
		}
            }
    };

    // At this point, s12, s23 are non-positive. 
    // If s13 is positive, make s12, s23 positive, unless s12 or s23 = 0.
    vector<Int4> diag = {1,1,1};
    if (S[0][2] > 0)
    {
        if (S[0][1] == 0) diag[0] = -1;
        else if (S[1][2] == 0) diag[2] = -1;
        else {diag[0] = -1; diag[2] = -1;}
    
    };
    NRMat<Int4> arr = diagonal (diag);
    g = mprod (arr, g);
    S = mprod (mprod (arr, S), arr);

    //print_mat (S);
    assert (
    (S[0][0] <= S[1][1] && S[1][1] <= S[2][2]) &&
    (2 * abs (S[0][1]) <= S[0][0]*EPS) &&
    (2 * abs (S[0][2]) <= S[0][0]*EPS) &&
    (2 * abs (S[1][2]) <= S[1][1]*EPS) &&
    (((S[0][1] > 0) && (S[0][2] > 0) && (S[0][2] > 0)) |
    ((S[0][1] <= 0) && (S[0][2] <= 0) && (S[0][2] <= 0))));
    
    assert (-2 * (S[0][1] + S[0][2] + S[1][2]) <= (S[0][0] + S[1][1])*EPS);

    Double d1 = dist (S, mprod (mprod (g, S_input),
                                    transpose (g)));
    NRMat<Double> zeros (ndim, ndim, 0.0);
    Double d2 = dist (S_input, zeros);
    assert (d1 <= d2);

    //return forward_as_tuple (g, S);
}

#endif /*PUT_MINKOWSKI_REDUCED_LATTICE_HH_*/
