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
# include <map>
# include <cmath>
# include <iostream>
# include <iomanip>
# include "../utility_lattice_reduction/mlist.hh"
# include "../utility_data_structure/nrutil_nr.hh"
# include "../utility_lattice_reduction/error_stable_bravais_2D.hh"
# include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
# include "error_stable_bravais_3D.hh"
# include "../ControlParam.hh"
using namespace std;

const vector <string> BravaisLatticeDetermination::listNames = {
"triclinic", "primitive monoclinic",
"base-centered monoclinic", "primitive orthorhombic",
"base-centered orthorhombic", "body-centered orthorhombic",
"face-centered orthorhombic", "primitive tetragonal",
"body-centered tetragonal", "rhombohedral", "hexagonal",
"primitive cubic", "body-centered cubic",
"face-centered cubic"};

Double rad2deg (const Double& rad)
{
    return 180.0 * rad / M_PI;
}

inline vector<Double> S2lattice (const NRMat<Double>& S)
{
    Double a,b,c,alpha,beta,gamma;
    a = sqrt (S[0][0]); b = sqrt (S[1][1]); c = sqrt (S[2][2]);
    alpha = acos (S[1][2] / (b * c)); alpha = rad2deg (alpha);
    beta = acos (S[0][2] / (a * c)); beta = rad2deg (beta);
    gamma = acos (S[0][1] / (a * b)); gamma = rad2deg (gamma);
    return  {a, b, c, alpha, beta, gamma};
}

static const vector<vector<Int4>> hF = {{1,1,0},{1,-1,0},{1,1,2}};
static const vector<vector<Int4>> hI = {{1,1,-1},{1,-1,0},{0,0,1}};
static const vector<vector<Int4>> hR = {{-1,0,1},{0,1,-1},{-1,-1,-1}};
static const vector<vector<Int4>> hB = {{1,1,0},{1,-1,0},{0,0,1}};
const NRMat<Int4> BravaisLatticeDetermination::ConstantType::h_F = put_matrix (hF);
const NRMat<Int4> BravaisLatticeDetermination::ConstantType::h_I = put_matrix (hI);
const NRMat<Int4> BravaisLatticeDetermination::ConstantType::h_R = put_matrix (hR);
const NRMat<Int4> BravaisLatticeDetermination::ConstantType::h_B = put_matrix (hB);

const vector<NRMat<Int4>>& BravaisLatticeDetermination::ConstantType::putCmP()
{
    static const vector<NRMat<Int4>> ans = {identity_matrix<Int4> (3), put_perm_mat3 (0,2,1), put_perm_mat3 (1,0,2)};
    return ans;
}

const vector<NRMat<Int4>>& BravaisLatticeDetermination::ConstantType::putCmC (bool doesPrudentSearch)
{
	static const vector<NRMat<Int4>> ans[2] = { MatrixList::put_CmC_list(false), MatrixList::put_CmC_list(true) };
	return (doesPrudentSearch?ans[1]:ans[0]);
}


const vector<NRMat<Int4>>& BravaisLatticeDetermination::ConstantType::putCoF()
{
    static const vector<NRMat<Int4>> cof = {h_F, mprod(h_F, put_perm_mat3 (0,2,1)), mprod (h_F, put_perm_mat3 (1,2,0))};
    return cof;
}

const vector<NRMat<Int4>>& BravaisLatticeDetermination::ConstantType::putCoI()
{
	static const vector<NRMat<Int4>> coi = {h_I, mprod (h_I, put_perm_mat3 (0,2,1)), mprod (h_I, put_perm_mat3 (1,2,0))};
    return coi;
}


const vector<NRMat<Int4>>& BravaisLatticeDetermination::ConstantType::putChR (bool doesPrudentSearch)
{
	static const vector<NRMat<Int4>> ans[2] = { MatrixList::put_ChR_list(false), MatrixList::put_ChR_list(true) };
	return (doesPrudentSearch?ans[1]:ans[0]);
}


void BravaisLatticeDetermination::InputType::set(const ControlParam& cData)
{
    Sobs = cData.putSobs();
    eps = cData.putResolution();
    doesPrudentSearch = cData.DoesPrudentSymmetrySearch();
    axisForBaseCenteredSymmetry = cData.put_ABCaxisString();
    axisForRhombohedralSymmetry = cData.put_RHaxisString();
}


void BravaisLatticeDetermination::InputType::display()
{
    cout << "----S input----" << endl;
	print_mat (put_Sobs());
	cout << "---input parameters---" << endl;
	cout << "eps = " << put_eps() << " doesPrudentSearch = " << put_doesPrudentSearch() << endl;
	cout << "axisBaseCentered = " << put_axisForBaseCenteredSymmetry() << endl;
	cout << "axisRhombohedral = " << put_axisForRhombohedralSymmetry() << endl;
}


void BravaisLatticeDetermination::reset()
{
    for (string name : listNames)
    {
        m_objResult[name].clear();
    }
}

void BravaisLatticeDetermination::to_output_format(vector<string>& namesOrder, map<string,vector<vector<Double>>>& arg) const
{
	namesOrder.clear();
	arg.clear();
    vector<Double> lattice_info;
    for (auto it=listNames.rbegin(); it<listNames.rend(); it++)
    {
        const string ctype = key2str(*it, m_objInput.put_axisForBaseCenteredSymmetry());
        namesOrder.push_back(ctype);
        const map<string, vector<ResultType>>::const_iterator itmap = m_objResult.find(*it);
        if( itmap == m_objResult.end() ) continue;

        vector< vector<Double> > lvec;
        for (auto it2=itmap->second.begin(); it2 < itmap->second.end(); it2++)
        {
            lattice_info = S2lattice (it2->S);
            lattice_info.push_back (it2->distance);
            lvec.push_back (lattice_info);
        }
        arg.insert (make_pair (ctype, lvec));
    }
}


vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::primitive_monoclinic(
                    const NRMat<Double>& S_obs, const NRMat<Int4>& gb, const Double& eps)
{
    vector<ResultType> ans;
    NRMat<Int4> h2; NRMat<Double> S2;
    for (NRMat<Int4> g : ConstantType::putCmP())
    {
        const NRMat<Int4> g2 = mprod (g, gb);
        const NRMat<Double> S_new = mprod (mprod (g2, S_obs), transpose (g2));

        const vector<vector<Double>> Sv = {
                {S_new[0][0],         0.0, S_new[0][2]},
                {        0.0, S_new[1][1],         0.0},
                {S_new[0][2],         0.0, S_new[2][2]}};

        const NRMat<Double> S = put_matrix (Sv);
        if (check_equiv (S_new, S, eps))
        {
        	const vector<vector<Double>> S2_vec
        	    = { {S[0][0], S[0][2]},
                    {S[2][0], S[2][2]}};
        	const NRMat<Double> S2_ = put_matrix (S2_vec);
            
            gauss_algorithm (S2_, h2, S2);
            
            const vector<vector<Int4>> h_vec
                  = {  {h2[0][0], 0, h2[0][1]},
                    {          0, 1,        0},
                    {   h2[1][0], 0, h2[1][1]} };
            const NRMat<Int4> h = put_matrix (h_vec);
            ResultType gs;
            gs.g = mprod (h, g2);
            gs.S = mprod (mprod (h, S), transpose (h));
            ans.push_back (gs);
        }
    }
    return ans;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::base_centered_monoclinic (
                    const NRMat<Double>& S_obs,
                    const NRMat<Int4>& gd, const Double& eps,
                    const bool& doesPrudentSearch,
                    const string& axisForBaseCenteredSymmetry)
{
    assert (axisForBaseCenteredSymmetry == "A" || axisForBaseCenteredSymmetry == "B" || axisForBaseCenteredSymmetry == "C");
    const vector<NRMat<Int4>>& CmCvec = ConstantType::putCmC(doesPrudentSearch);

    vector<ResultType> ans;
    NRMat<Int4> h, h2;
    NRMat<Double> S2;
    for (NRMat<Int4> g : CmCvec)
    {
    	const NRMat<Int4> g2 = mprod (g, gd);
    	const NRMat<Double> S_new = mprod (mprod (g2, S_obs), transpose (g2));
    	const vector<vector<Double>> S_vec
    	      = {  {S_new[0][0],           0, S_new[0][2]},
                {          0, S_new[1][1],           0},
                {   S_new[0][2],           0, S_new[2][2]} };
    	const NRMat<Double> S = put_matrix (S_vec);
        if (check_equiv (S_new, S, eps))
        {
        	const vector<vector<Double>> S2_vec
        	      = { {S[0][0],S[0][2]},
                      {S[2][0],S[2][2]}};
            const NRMat<Double> S2_ = put_matrix (S2_vec);
            gauss_algorithm (S2_, h2, S2);
            if (h2[1][0] % 2 != 0)
            {
                if (h2[0][0] % 2 != 0)
                {
                	static const NRMat<Int4> arrv = put_matrix(vector<vector<Int4>>({{1,0},{-1,-1}}));
                	h2 = mprod (arrv, h2);
                }
                else
                {
                	static const NRMat<Int4> arrv = put_matrix(vector<vector<Int4>>({{0,1},{1,0}}));
                	h2 = mprod (arrv, h2);
                }
            }
            const vector<vector<Int4>> h_vec
                  = {  {h2[0][0], 0, h2[0][1]},
                    {       0, 1,        0},
                    {h2[1][0], 0, h2[1][1]}};
            h = put_matrix (h_vec);
            if (axisForBaseCenteredSymmetry == "A")
            {
            	h = mprod (put_perm_mat3 (1,2,0), h);
            }
            else if (axisForBaseCenteredSymmetry == "C")
            {
            	h = mprod (put_perm_mat3 (2,0,1), h);
            }

            ResultType gs;
            gs.g = mprod (h, g2);
            gs.S = mprod (mprod (h, S), transpose (h));
            ans.push_back (gs);
        }
    }
    return ans;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::face_centered_orthorhombic (
                        const NRMat<Double>& S_obs, const NRMat<Int4>& gd,
                        const Double& eps)
{
    vector<ResultType> ans;
    for (NRMat<Int4> g : ConstantType::putCoF())
    {
    	const NRMat<Int4> g2 = mprod (g, gd);
    	const NRMat<Double> S_new = mprod (mprod (g2, S_obs), transpose (g2));
    	const vector<vector<Double>> S_vec
    	   = {  {S_new[0][0],0,0},
                {0,S_new[1][1],0},
                {0,0,S_new[2][2]}};
        const NRMat<Double> S = put_matrix (S_vec);
        
        if (check_equiv (S_new, S, eps))
        {
        	const vector<Double> diagS = {S[0][0], S[1][1], S[2][2]};
            const vector<Int4> indices = argsort (diagS);
            const NRMat<Int4> gp = put_perm_mat3 (indices[0], indices[1], indices[2]);
            ResultType gs;
            gs.g = mprod (gp, g2);
            gs.S = mprod (mprod (gp, S), transpose (gp));
            ans.push_back (gs);
        }
    }
    return ans;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::body_centered_orthorhombic (
    const NRMat<Double>& S_obs, const NRMat<Int4>& gd2, const Double& eps)
{
    vector<ResultType> ans;
    vector<Int4> indices;
    NRMat<Int4> g2, gp, g2_; NRMat<Double> S, S_new;
    vector<vector<Double>> Sv; vector<Double> diagS;
    
    for (NRMat<Int4> g : ConstantType::putCoI())
    {
        g2 = mprod (g, gd2);
        S_new = mprod (mprod (g2, S_obs), transpose (g2));
        Sv = {  {S_new[0][0],0,0},
                {0,S_new[1][1],0},
                {0,0,S_new[2][2]}};
        S = put_matrix (Sv);

        if (check_equiv (S_new, S, eps))
        {
            diagS = {S[0][0], S[1][1], S[2][2]};
            indices = argsort (diagS);
            //gp = permute_matrix (indices);
            gp = put_perm_mat3 (indices[0], indices[1], indices[2]);
            ResultType gs;
            gs.g = mprod(gp, g2);
            gs.S = mprod (mprod (gp, S), transpose (gp));
            ans.push_back (gs);
        }
    }
    return ans;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::rhombohedral_centering (
    const NRMat<Double>& S_obs, const NRMat<Int4>& gd, const Double& eps,
    const bool& doesPrudentSearch, const string& axisForRhombohedralSymmetry)
{
    assert (axisForRhombohedralSymmetry == "Rhombohedral" || axisForRhombohedralSymmetry == "Hexagonal");
    static const NRMat<Int4>& h_R = ConstantType::put_hR();
    const vector<NRMat<Int4>>& ChRvec = ConstantType::putChR(doesPrudentSearch);

    vector<ResultType> ans;
    for (NRMat<Int4> g : ChRvec)
    {
        const NRMat<Int4> g2 = mprod (g, gd);
        const NRMat<Double> S_new = mprod (mprod (g2, S_obs), transpose (g2));

        const Double a = (S_new[0][0] + S_new[1][1] + S_new[2][2]) / 3.;
        const Double d = (S_new[0][1] + S_new[0][2] + S_new[1][2]) / 3.;
        const vector<vector<Double>> S_vec = {{a,d,d},{d,a,d},{d,d,a}};
        const NRMat<Double> S = put_matrix (S_vec);

        if (check_equiv (S_new, S, eps))
        {
            ResultType gs;
            if (axisForRhombohedralSymmetry == "Hexagonal")
            {
                gs.g = mprod (h_R, g2);
                gs.S = mprod (mprod (h_R, S), transpose (h_R));
            }
            else
            {
            	gs.g = g2;
            	gs.S = S;
            }
            ans.push_back (gs);
        }
    }
    return ans;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::primitive_monoclinic_to_orthorhombic (
                                  const vector<BravaisLatticeDetermination::ResultType>& mlist, const Double& eps)
{
    vector<ResultType> ans;
    for (ResultType gsm : mlist)
    {
    	const NRMat<Int4>& g = gsm.g;
    	const NRMat<Double>& Sm = gsm.S;
        const Double& s11 = Sm[0][0];
        const Double& s22 = Sm[1][1];
        const Double& s33 = Sm[2][2];
        const vector<vector<Double>> So_vec = {{s11,0,0},{0,s22,0},{0,0,s33}};
        const NRMat<Double> So = put_matrix (So_vec);
        if (check_equiv (Sm, So, eps))
        {
        	const vector<Double> diagSo = {So[0][0], So[1][1], So[2][2]};
        	const vector<Int4> indices = argsort (diagSo);
        	const NRMat<Int4> gp = put_perm_mat3 (indices[0], indices[1], indices[2]);
            ResultType gs;
            gs.g = mprod (gp, g);
            gs.S = mprod (mprod (gp, So), transpose (gp));
            ans.push_back (gs);
        }
    }
    return ans;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::base_monoclinic_to_orthorhombic (
                        const vector<BravaisLatticeDetermination::ResultType>& mlist, const Double& eps,
                        const string& axisForBaseCenteredSymmetry)
{
    NRMat<Int4> gp, g_new; NRMat<Double> So;
    vector<ResultType> ans;
    for (ResultType gsm : mlist)
    {
        const NRMat<Int4>& g = gsm.g;
        const NRMat<Double>& Sm = gsm.S;
        const Double& s11 = Sm[0][0];
        const Double& s22 = Sm[1][1];
        const Double& s33 = Sm[2][2];
        const vector<vector<Double>> So_vec = {{s11,0,0}, {0,s22,0}, {0,0,s33}};
        So = put_matrix (So_vec);
        if (check_equiv (Sm, So, eps))
        {
            if (axisForBaseCenteredSymmetry == "A" || axisForBaseCenteredSymmetry == "C")
            {
                if (axisForBaseCenteredSymmetry == "A")
                {
                	gp = put_perm_mat3 (2,0,1);
                }
                else
                {
                	gp = put_perm_mat3 (1,2,0);
                }
                g_new = mprod (gp, g);
                So = mprod (mprod (gp, So), transpose (gp));
            }
            else g_new = g;
            
            ResultType gs;
            if (So[0][0] > So[1][1])
            {
                static const NRMat<Int4> gp =  put_perm_mat3 (1,0,2);
                gs.g = mprod (gp, g_new);
                gs.S = mprod (mprod (gp, So), transpose (gp));
            }
            else
            {
            	gs.g = g_new;
            	gs.S = So;
            }
            ans.push_back (gs);
        }
    }
    return ans;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::orthorhombic_to_tetragonal (
                           const vector<BravaisLatticeDetermination::ResultType>& mlist, const Double& eps)
{
    static const vector<vector<Int4>> ijks = {{0,1,2},{0,2,1},{1,2,0}};

    vector<ResultType> ans;
    for (ResultType gsm : mlist)
    {
    	const NRMat<Int4>&  g = gsm.g;
    	const NRMat<Double>& So = gsm.S;
        for (vector<Int4> ijk : ijks)
        {
            const Int4 i = ijk[0];
            const Int4 j = ijk[1];
            const Int4 k = ijk[2];
            const Double s = (So[i][i] + So[j][j]) / 2.0;
            const vector<vector<Double>> St_vec = {{s,0,0},{0,s,0},{0,0,So[k][k]}};
            const NRMat<Double> St = put_matrix (St_vec);
            const NRMat<Int4> gp = put_perm_mat3 (i,j,k);
            const NRMat<Int4> g_new = mprod(gp, g);
            if (check_equiv (mprod (mprod (gp, So), transpose (gp)), St, eps))
            {
                ResultType gs;
                gs.g = g_new;
                gs.S = St;
                ans.push_back (gs);
            }
        }
    }
    return ans;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::primitive_monoclinic_to_hexagonal (
                                                const vector<ResultType>& ans_mP, const Double& eps)
{
    static const NRMat<Int4> gp = put_perm_mat3 (0,2,1);
    vector<ResultType> ans_hP;
    for (ResultType gsm : ans_mP)
    {
        const NRMat<Int4>& g = gsm.g;
        const NRMat<Double>& S_mP = gsm.S;
        const Double& s11 = S_mP[0][0];
        const Double& s22 = S_mP[1][1];
        const Double& s33 = S_mP[2][2];
        const Double& s13 = S_mP[0][2];
        const Double s = s11 + s33 + 2 * s13;
        const vector<vector<Double>> S_hP_vec = {{s, -s/2, 0}, {-s/2, s, 0}, {0, 0, s22}};
        const NRMat<Double> S_hP = put_matrix (S_hP_vec);
        const NRMat<Int4> g_new = mprod (gp, g);
        if( check_equiv (mprod (mprod (gp, S_mP), transpose (gp)), S_hP, eps) )
        {
        	ResultType gs;
        	gs.g = g_new;
        	gs.S = S_hP;
        	ans_hP.push_back (gs);
        }
    }
    return ans_hP;
}

vector<BravaisLatticeDetermination::ResultType> BravaisLatticeDetermination::orthorhombic_to_cubic (
                                            const vector<ResultType>& mlist, const Double& eps)
{
    vector<ResultType> ans;
    for (ResultType gsm : mlist)
    {
        const NRMat<Int4>& g = gsm.g;
        const NRMat<Double>& So = gsm.S;
        const Double s = (So[0][0] + So[1][1] + So[2][2]) / 3.0;
        const vector<vector<Double>> Sc_vec = {{s, 0, 0}, {0, s, 0}, {0, 0, s}};
        const NRMat<Double> Sc = put_matrix (Sc_vec);
        if (check_equiv (So, Sc, eps))
        {
        	ResultType gs;
        	gs.g = g;
        	gs.S = Sc;
            ans.push_back (gs);
        }
    }
    return ans;
}

void BravaisLatticeDetermination::set_bravais_class (const InputType& input)
{
  m_objInput = input;
  this->reset();
  const NRMat<Double> Sobs = m_objInput.put_Sobs();
  const Double eps = m_objInput.put_eps();
  const bool doesPrudentSearch = m_objInput.put_doesPrudentSearch();
  const string axisForBasecenteredSymmetry = m_objInput.put_axisForBaseCenteredSymmetry();
  const string axisForRhombohedralSymmetry = m_objInput.put_axisForRhombohedralSymmetry();

  NRMat<Double> S_buer, S_del, S_del_inv;
  NRMat<Int4> gb, gd, gd2;
  buerger_reduction (Sobs, gb, S_buer);
  Delaunay_reduction (S_buer, gd, S_del); gd = mprod (gd, gb);
  Delaunay_reduction_of_inverse (S_buer, gd2, S_del_inv); gd2 = mprod (gd2, gb);

  ResultType gstri;
  gstri.g = gb;
  gstri.S = S_buer;
  m_objResult["triclinic"].resize(1, gstri);
  m_objResult["primitive monoclinic"] = primitive_monoclinic (Sobs, gb, eps);
  m_objResult["base-centered monoclinic"] = base_centered_monoclinic (Sobs, gd, eps, doesPrudentSearch,
                                          	  	  	  	  	  	  	  	  	axisForBasecenteredSymmetry);
  m_objResult["face-centered orthorhombic"] = face_centered_orthorhombic (Sobs, gd, eps);
  m_objResult["body-centered orthorhombic"] = body_centered_orthorhombic (Sobs, gd2, eps);
  m_objResult["rhombohedral"] = rhombohedral_centering (Sobs, gd, eps, doesPrudentSearch,
                                        									axisForRhombohedralSymmetry);
  m_objResult["primitive orthorhombic"] = primitive_monoclinic_to_orthorhombic (
                								m_objResult["primitive monoclinic"], eps );
  m_objResult["base-centered orthorhombic"] = base_monoclinic_to_orthorhombic (
                								m_objResult["base-centered monoclinic"], eps, axisForBasecenteredSymmetry );
  m_objResult["primitive tetragonal"] = orthorhombic_to_tetragonal (
                								m_objResult["primitive orthorhombic"], eps);

  // Body-centered orthorombic --> Body-centered tetragonal
  m_objResult["body-centered tetragonal"] = orthorhombic_to_tetragonal (
                    							m_objResult["body-centered orthorhombic"], eps);
  // Primitive Monoclinic --> Hexagonal
  m_objResult["hexagonal"] = primitive_monoclinic_to_hexagonal (
                    							m_objResult["primitive monoclinic"], eps);

  // Primtive orthorhombic --> Primitive cubic
  m_objResult["primitive cubic"] = orthorhombic_to_cubic (
                    							m_objResult["primitive orthorhombic"], eps);

  // Body-centered orthorhombic --> Body-centered cubic
  m_objResult["body-centered cubic"] = orthorhombic_to_cubic (
                    							m_objResult["body-centered orthorhombic"], eps);

  // Face-centered orthorhombic -->  Face-centered cubic
  m_objResult["face-centered cubic"] = orthorhombic_to_cubic (
		  	  	  	  	  	  	  	  	  	  	m_objResult["face-centered orthorhombic"], eps);

  // Set distances and sort.
  for(auto it=m_objResult.begin(); it!=m_objResult.end(); it++)
  {
	for(auto it2=it->second.begin(); it2<it->second.end(); it2++)
	{
        it2->setDistance(Sobs);
	}
    sort(it->second.begin(), it->second.end());
  }
}


static string print_matrices(const NRMat<Int4>& g, const NRMat<Double>& S,
                const NRMat<Double>& S2, const Double& d)
{
    string text = "--- g, projection of g*Sobs*g^T, g*Sobs*g^T ---";
    text += "distance = " + to_string (d) + "\n";
    const Int4 nrow = g.nrows(), ncol = g.ncols();
    for (Int4 i = 0; i < nrow; i++)
    {
        text += "| ";
        for (Int4 j = 0; j < ncol; j++)
        {
            text += to_string (g[i][j]);
            if (j < ncol - 1) text += " , ";
            else text += " | ";
        }

        for (Int4 j = 0; j < ncol; j++)
        {
            text += to_string (S[i][j]);
            if (j < ncol - 1) text += " , ";
            else text += " | ";
        }

        for (Int4 j = 0; j < ncol; j++)
        {
            text += to_string (S2[i][j]);
            if (j < ncol - 1) text += " , ";
            else if (i == nrow - 1) text += " |";
            else text += " |\n";
        }
    }
    return text;
}

string BravaisLatticeDetermination::toText0() const
{
    NRMat<Int4> g;
    NRMat<Double> S, Sproj;
    stringstream strs;
    for (string bname : listNames)
    {
        const string bname_str = this->key2str (bname, this->m_objInput.put_axisForBaseCenteredSymmetry());
        const map<string, vector<ResultType>>::const_iterator it = m_objResult.find(bname);
        if( it == m_objResult.end() ) continue;

        strs <<"---- " + bname_str +  "  length : " + to_string (it->second.size()) + " ----" << endl;

        for (auto it2=it->second.begin(); it2 < it->second.end(); it2++)
        {
            const NRMat<Int4>& g = it2->g;
            const NRMat<Double>& Sproj = it2->S;
            S = mprod (mprod (g, m_objInput.put_Sobs()), transpose (g));
            strs << print_matrices (g, Sproj, S, it2->distance) << "\n";

        }
        strs << "" << endl;
    }
    return strs.str();
}

inline string float2string (const Double& x)
{
    stringstream strs;
    strs << scientific << setprecision (4) << x;
    return strs.str();
}

inline string adjust_name_string (const string& name, const Int4& maxLen)
{
    stringstream strs;
    strs << right << setw(maxLen) << fixed << name;
    return strs.str();
}

static string make_text_summary (const vector<string>& namesOrder, const map<string, vector<vector<Double>>>& bmap)
{
    string text = "   <!-- << Unit-cell Parameters with the Minimal Distance >>\n";
    text += "   Bravais_type  Number_of_candidates, a, b, c, alpha, beta, gamma, dstance_from_unit_cell\n";

    static const Int4 maxLen = 15;
    for (string ctype : namesOrder)
    {
    	map<string, vector<vector<Double>>>::const_iterator it=bmap.find(ctype);
    	if( it == bmap.end() ) continue;
        const vector<vector<Double>>& lattices_info = it->second;

        ctype = adjust_name_string (ctype, maxLen);
        text += "   " + ctype + "   ";

        if( !(lattices_info.empty()) )
        {
            text += to_string(lattices_info.size()) + "   ";
            for (Double p : lattices_info[0])
            {
            	text += "   "s + float2string (p);
            }
        }
        else text += "0";
		text += "\n";
    }
    text += "   -->\n";
    return text;
}

static string make_text_all_candidates (
                const vector<vector<Double>>& lattices_info,
                string ctype)
{
    string text = "  <!-- Candidates for " + ctype + " -->" + "\n";
    if ( lattices_info.empty() ) return text;

    text += "  <UnitCellCandidate>\n";
        text += "   <CrystalSystem> " + ctype + " </CrystalSystem>\n";
        text += "   <UnitCellParameters>\n";
        text += "   <!-- a, b, c, alpha, beta, gamma, distance from input cell -->\n";
        for (auto it = lattices_info.begin(); it < lattices_info.end(); it++)
        {
            text += "    ";
            for (Double l : *it) text += "   " + float2string (l);
            text += "\n";
        }
        text += "   </UnitCellParameters>\n";
    text += "  </UnitCellCandidate>\n";
    return text;
}

string BravaisLatticeDetermination::toText() const
{
    map<string, vector<vector<Double>>> bmap;
    vector<string> namesOrder;
    this->to_output_format(namesOrder, bmap);

    string text = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    text += "<ZCodeParameters>\n";
    text += " <BLDConographOutput>\n";
    text += make_text_summary (namesOrder, bmap);
    for (string ctype : namesOrder)
    {
    	auto it = bmap.find(ctype);
    	if( it == bmap.end() ) continue;
        text += make_text_all_candidates (it->second, ctype);
        text += "\n";
    }

    text += " </BLDConographOutput>\n";
    text += "</ZcodeParameters>";

    return text;
}
