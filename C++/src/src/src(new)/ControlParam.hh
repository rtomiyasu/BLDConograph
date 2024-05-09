/*
 * The MIT License

   BLDConograph (Bravais lattice determination module in Conograph)

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
#ifndef _ControlParam_h_
#define _ControlParam_h_
// ControlParam.hh

# include <iostream>
#include <map>
#include "RietveldAnalysisTypes.hh"
#include "bravais_type/enumAxis.hh"
#include "zerror_type/error_out.hh"
#include "utility_data_structure/VecDat3.hh"
#include "utility_data_structure/SymMat.hh"
#include "utility_func/lattice_constant.hh"
#include "utility_rw_param/I_ReadData.hh"
using namespace std;

class ControlParam : public I_ReadData
{
private:
	enum{ NUM_LS = 14 };
	static const map<eABCaxis, string> ABCaxisLabel;
	static const map<eRHaxis, string> RHaxisLabel;
	
	// Parameters for search.
	static const pair< RWParamProperty, RWParamData<void*> > LatticeParameter_Data;
	static const pair< RWParamProperty, RWParamData<Double> > LatticeParameters_Data[6];

	static const pair< RWParamProperty, RWParamData<bool> > DoesPrudentSymSearch_Data;
	static const pair< RWParamProperty, RWParamData<bool> > OutputSymmetry_Data[NUM_LS];

	static const pair< RWParamProperty, RWParamData<Double> > Resol_Data;	// unit->Angstrom.
		
//	static const pair< RWParamProperty, RWParamData<Int4> > NumCores_Data;
	static const pair< RWParamProperty, RWParamData<string> > MonoBaseAxis_Data;
	static const pair< RWParamProperty, RWParamData<string> > RhomAxis_Data;

	VecDat3<Double> m_length;
	VecDat3<Double> m_angle;

	bool DoesPrudentSymSearch;
	bool OutputSymmetry[NUM_LS];

	Double Resol;	// unit->Angstrom.
	Double eps;

	// Enviromental parameters.
	static const Int4 NumCores;
	string MonoBaseAxis;
	string RhomAxis;

protected:
	ZErrorMessage checkData(const RWParam_void& param) const;

public:
    ControlParam();
    virtual ~ControlParam();

	// Set functions.
    // Parameters for search.
    inline void setDoesPrudentSeartch(const bool& arg){ DoesPrudentSymSearch = arg; };
	inline void setOutputSymmetry(const bool* arg){ for(Int4 i=0; i<NUM_LS; i++) OutputSymmetry[i] = arg[i]; };

	inline void setResolution(const Double& arg) { Resol = arg; };

	// Put functions.
	inline const bool& DoesPrudentSymmetrySearch() const { return DoesPrudentSymSearch; };
//	inline const bool& putOutputSymmetry(const eBravaisType& i) const { return OutputSymmetry[(Int4)i]; };

	inline const Double& putResolution() const { return Resol; };

	inline const Int4& putNumberOfThreadsToUse() const { return NumCores; };
	inline eRHaxis putRhombohedralAxis() const { return find_key(RHaxisLabel, RhomAxis); };
	inline eABCaxis putBaseCenteredAxis() const { return find_key(ABCaxisLabel, MonoBaseAxis); };

	const string& putTagLabel() const;
    void setData(const RWParamProperty& parent_prop,
			vector<RWParam_void>& tray);
	inline const string& put_ABCaxisString () const {return MonoBaseAxis;};
	inline const string& put_RHaxisString () const {return RhomAxis;};
	inline NRMat<Double> putSobs () const;
};

template<class T>
inline NRMat<T> symMat2nrMat (SymMat<T> m)
{
	Int4 ndim = m.size();
	NRMat<T> ans (ndim, ndim, 0);
	for (Int4 i = 0; i < ndim; i++)
	{
		ans[i][i] = m(i,i);
		for (Int4 j = i + 1; j < ndim; j++)
		{
			ans[i][j] = ans[j][i] = m(j,i);
		}
	}
	return ans;
}

inline NRMat<Double> ControlParam::putSobs () const
{
	SymMat<Double> Sobs_inv (3), S_super (4);
	NRMat<Double> Sobs;
	NRMat<Int4> trans_mat(4,3);
	calCoParameter(m_length, m_angle, Sobs_inv);
	Sobs = inverse_mat_3x3 (symMat2nrMat (Sobs_inv));
	return Sobs;
}

#endif
