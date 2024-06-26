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

#ifndef MLIST_HH_
#define MLIST_HH_
#include <iostream>
#include <cassert>
#include <stdio.h>
#include "../utility_data_structure/nrutil_nr.hh"

using namespace std;

class MatrixList
{
private:
    enum{ lenCmC_1 = 21, lenCmC = 69 };
    enum{ lenChR_1 = 16, lenChR = 64 };

    static const vector<NRMat<Int4>> CmC_list;
    static const vector<NRMat<Int4>> ChR_list;

public:
    static vector<NRMat<Int4>> put_CmC_list (bool doesPrudentSearch);
    static vector<NRMat<Int4>> put_ChR_list (bool doesPrudentSearch);
};


#endif

/*Int4 len_CmC_1 = 21; Int4 len_CmC = 69;
Int4 len_ChR_1 = 16; Int4 len_ChR = 64;

inline vector<vector<vector<Int4>>> cmcList = {
    // CmC_1 (21 arrays)
    {{1, 1, 0}, {1, -1, 0}, {0, 0, 1}},
    {{1, 0, 1}, {1, 0, -1}, {0, 1, 0}},
    {{1, 0, -1}, {1, 0, 1}, {0, 1, 1}},
    {{1, -1, 0}, {1, 1, 0}, {0, 1, 1}},
    {{2, 1, 0}, {0, -1, 0}, {0, 0, 1}},
    {{2, 1, 0}, {0, -1, 0}, {-1, -1, -1}},
    {{2, 0, 1}, {0, 0, -1}, {0, 1, 0}},
    {{2, 0, 1}, {0, 0, -1}, {-1, -1, -1}}, //8
    {{1, -1, -1}, {1, 1, 1}, {0, 1, 0}},
    {{1, -1, -1}, {1, 1, 1}, {0, 0, 1}},
    {{2, 1, 1}, {0, -1, -1}, {-1, 0, -1}},
    {{0, -1, -1}, {2, 1, 1}, {0, 1, 0}},
    {{0, 1, 1}, {0, 1, -1}, {1, 0, 0}},
    {{1, 2, 0}, {-1, 0, 0}, {0, 0, 1}},
    {{1, 2, 0}, {-1, 0, 0}, {-1, -1, -1}},
    {{0, 2, 1}, {0, 0, -1}, {1, 0, 0}}, //16
    {{-1, 1, -1}, {1, 1, 1}, {1, 0, 0}},
    {{-1, 0, -1}, {1, 2, 1}, {1, 0, 0}},
    {{1, 0, 2}, {-1, 0, 0}, {0, 1, 0}},
    {{0, 1, 2}, {0, -1, 0}, {1, 0, 0}},
    {{-1, -1, 0}, {1, 1, 2}, {1, 0, 0}},

    // CmC_2 (48 arrays)
    {{1, 1, 0}, {1, -1, 0}, {1, 0, 1}},
    {{1, 1, 0}, {1, -1, 0}, {0, -1, -1}},
    {{1, 0, 1}, {1, 0, -1}, {1, 1, 0}}, //24
    {{1, 0, 1}, {1, 0, -1}, {0, -1, -1}},
    {{1, 0, -1}, {1, 0, 1}, {0, 1, 0}},
    {{1, 0, -1}, {1, 0, 1}, {-1, -1, -1}},
    {{1, -1, 0}, {1, 1, 0}, {0, 0, 1}},
    {{1, -1, 0}, {1, 1, 0}, {-1, -1, -1}},
    {{1, 1, 1}, {1, -1, -1}, {0, 0, -1}},
    {{1, 1, 1}, {1, -1, -1}, {0, -1, 0}},
    {{0, 0, -1}, {2, 0, 1}, {0, -1, 0}}, //32
    {{0, 0, -1}, {2, 0, 1}, {1, 1, 1}},
    {{0, -1, 0}, {2, 1, 0}, {0, 0, -1}},
    {{0, -1, 0}, {2, 1, 0}, {1, 1, 1}},
    {{2, 1, 1}, {0, -1, -1}, {0, 1, 0}},
    {{2, 1, 1}, {0, -1, -1}, {0, 0, 1}},
    {{0, -1, -1}, {2, 1, 1}, {1, 1, 0}},
    {{0, -1, -1}, {2, 1, 1}, {1, 0, 1}},
    {{0, 1, 1}, {0, 1, -1}, {1, 1, 0}}, //40
    {{0, 1, 1}, {0, 1, -1}, {-1, 0, -1}},
    {{0, 1, -1}, {0, 1, 1}, {1, 0, 0}},
    {{0, 1, -1}, {0, 1, 1}, {-1, -1, -1}},
    {{1, 1, 1}, {-1, 1, -1}, {0, 0, -1}},
    {{1, 1, 1}, {-1, 1, -1}, {-1, 0, 0}},
    {{0, 0, -1}, {0, 2, 1}, {-1, 0, 0}},
    {{0, 0, -1}, {0, 2, 1}, {1, 1, 1}},
    {{-1, 0, 0}, {1, 2, 0}, {0, 0, -1}}, //48
    {{-1, 0, 0}, {1, 2, 0}, {1, 1, 1}},
    {{1, 2, 1}, {-1, 0, -1}, {1, 0, 0}},
    {{1, 2, 1}, {-1, 0, -1}, {0, 0, 1}},
    {{-1, 0, -1}, {1, 2, 1}, {1, 1, 0}},
    {{-1, 0, -1}, {1, 2, 1}, {0, 1, 1}},
    {{1, 1, 1}, {-1, -1, 1}, {0, -1, 0}},
    {{1, 1, 1}, {-1, -1, 1}, {-1, 0, 0}},
    {{0, -1, 0}, {0, 1, 2}, {-1, 0, 0}}, //56
    {{0, -1, 0}, {0, 1, 2}, {1, 1, 1}},
    {{-1, 0, 0}, {1, 0, 2}, {0, -1, 0}},
    {{-1, 0, 0}, {1, 0, 2}, {1, 1, 1}},
    {{1, 1, 2}, {-1, -1, 0}, {1, 0, 0}},
    {{1, 1, 2}, {-1, -1, 0}, {0, 1, 0}},
    {{-1, -1, 0}, {1, 1, 2}, {1, 0, 1}},
    {{-1, -1, 0}, {1, 1, 2}, {0, 1, 1}},
    {{1, 0, 0}, {1, 2, 2}, {0, 1, 0}}, //64
    {{1, 0, 0}, {1, 2, 2}, {0, 0, 1}},
    {{0, 1, 0}, {2, 1, 2}, {1, 0, 0}},
    {{0, 1, 0}, {2, 1, 2}, {0, 0, 1}},
    {{0, 0, 1}, {2, 2, 1}, {1, 0, 0}},
    {{0, 0, 1}, {2, 2, 1}, {0, 1, 0}} //69
                                            };

inline vector<vector<vector<Int4>>> chrList = {

    // ChR_1 (16 arrays)
    {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
    {{1, 0, 0}, {0, 1, 0}, {-1, -1, -1}},
    {{1, 0, 0}, {0, 0, 1}, {-1, -1, -1}},
    {{1, 0, 0}, {0, 0, -1}, {1, 1, 0}},
    {{1, 0, 0}, {0, 0, -1}, {0, -1, -1}},
    {{1, 0, 0}, {0, -1, 0}, {1, 0, 1}},
    {{1, 0, 0}, {0, -1, 0}, {0, -1, -1}},
    {{1, 0, 0}, {1, 1, 1}, {1, 1, 0}}, //8
    {{1, 0, 0}, {1, 1, 1}, {1, 0, 1}},
    {{0, 1, 0}, {0, 0, 1}, {-1, -1, -1}},
    {{0, 1, 0}, {0, 0, -1}, {1, 1, 0}},
    {{0, 1, 0}, {0, 0, -1}, {-1, 0, -1}},
    {{0, 1, 0}, {1, 1, 1}, {1, 1, 0}},
    {{0, 1, 0}, {1, 1, 1}, {0, 1, 1}},
    {{0, 0, 1}, {1, 1, 1}, {1, 0, 1}},
    {{0, 0, 1}, {1, 1, 1}, {0, 1, 1}}, // 16 (end of ChR_1)

    // ChR (48 arrays)
    {{1, 0, 0}, {0, 0, -1}, {0, 1, 1}},
    {{1, 0, 0}, {0, 0, -1}, {-1, -1, 0}},
    {{1, 0, 0}, {0, -1, 0}, {0, 1, 1}},
    {{1, 0, 0}, {0, -1, 0}, {-1, 0, -1}},
    {{1, 0, 0}, {-1, 0, -1}, {1, 1, 1}},
    {{1, 0, 0}, {-1, -1, 0}, {1, 1, 1}},
    {{0, 1, 0}, {0, 0, -1}, {1, 0, 1}},
    {{0, 1, 0}, {0, 0, -1}, {-1, -1, 0}}, //24
    {{0, 1, 0}, {0, -1, -1}, {1, 1, 1}},
    {{0, 1, 0}, {-1, -1, 0}, {1, 1, 1}},
    {{0, 0, 1}, {0, -1, -1}, {1, 1, 1}},
    {{0, 0, 1}, {-1, 0, -1}, {1, 1, 1}},
    {{1, 0, 0}, {0, 1, 0}, {0, 0, -1}},
    {{1, 0, 0}, {0, 1, 0}, {1, 0, 1}},
    {{1, 0, 0}, {0, 1, 0}, {0, -1, -1}},
    {{1, 0, 0}, {0, 1, 0}, {1, 1, 1}}, //32
    {{1, 0, 0}, {0, 0, 1}, {0, -1, 0}},
    {{1, 0, 0}, {0, 0, 1}, {1, 1, 0}},
    {{1, 0, 0}, {0, 0, 1}, {0, -1, -1}},
    {{1, 0, 0}, {0, 0, 1}, {1, 1, 1}},
    {{1, 0, 0}, {0, 0, -1}, {-1, -1, -1}},
    {{1, 0, 0}, {0, -1, 0}, {-1, -1, -1}},
    {{1, 0, 0}, {1, 1, 0}, {-1, -1, -1}},
    {{1, 0, 0}, {1, 0, 1}, {-1, -1, -1}}, //40
    {{1, 0, 0}, {0, 1, 0}, {0, 1, 1}},
    {{1, 0, 0}, {0, 1, 0}, {-1, 0, -1}},
    {{1, 0, 0}, {0, 0, -1}, {0, -1, 0}},
    {{0, 1, 0}, {0, 0, 1}, {1, 1, 0}},
    {{0, 1, 0}, {0, 0, 1}, {-1, 0, -1}},
    {{0, 1, 0}, {0, 0, 1}, {1, 1, 1}},
    {{0, 1, 0}, {0, 0, -1}, {-1, -1, -1}},
    {{1, 0, 0}, {0, -1, 0}, {1, 1, 1}}, //48
    {{0, 1, 0}, {1, 1, 0}, {-1, -1, -1}},
    {{0, 1, 0}, {0, 1, 1}, {-1, -1, -1}},
    {{1, 0, 0}, {0, 0, 1}, {0, 1, 1}},
    {{1, 0, 0}, {0, 0, 1}, {-1, -1, 0}},
    {{0, 1, 0}, {0, 0, 1}, {1, 0, 1}},
    {{0, 1, 0}, {0, 0, 1}, {-1, -1, 0}},
    {{0, 1, 0}, {0, 0, -1}, {1, 1, 1}},
    {{1, 0, 0}, {0, 0, -1}, {1, 1, 1}}, //56
    {{0, 0, 1}, {1, 0, 1}, {-1, -1, -1}},
    {{0, 0, 1}, {0, 1, 1}, {-1, -1, -1}},
    {{0, 0, 1}, {-1, 0, -1}, {-1, -1, -1}},
    {{0, 0, 1}, {0, -1, -1}, {-1, -1, -1}},
    {{0, 1, 0}, {-1, -1, 0}, {-1, -1, -1}},
    {{0, 1, 0}, {0, -1, -1}, {-1, -1, -1}},
    {{1, 0, 0}, {-1, -1, 0}, {-1, -1, -1}},
    {{1, 0, 0}, {-1, 0, -1}, {-1, -1, -1}} // 64
                                            };

inline vector<NRMat<Int4>> make_matrix_list (
                vector<vector<vector<Int4>>> matList,
                Int4 num)
{
    vector<NRMat<Int4>> ans_mats;

    if (num == 0) num = matList.size();

    for (int k = 0; k < num; k++)
    {
        ans_mats.push_back (put_matrix (matList[k]));
    };
    return ans_mats;
}

inline vector <NRMat<Int4>> CmC = make_matrix_list (cmcList, 0);
inline vector <NRMat<Int4>> ChR = make_matrix_list (chrList, 0);
inline vector <NRMat<Int4>> CmC_1 = make_matrix_list (cmcList, len_CmC_1);
inline vector <NRMat<Int4>> ChR_1 = make_matrix_list (chrList, len_ChR_1);


//vector<NRMat<Int4>> CmC = make_matrix_list (cmcList);
//vector<NRMat<Int4>> ChR = make_matrix_list (chrList);

//assert (CmC.size() == len_CmC);
//assert (ChR.size() == len_ChR);*/



