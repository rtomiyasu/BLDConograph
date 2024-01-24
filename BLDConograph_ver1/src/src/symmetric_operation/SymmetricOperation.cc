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
#include "SymmetricOperation.hh"
#include "../point_group/coset_representative_data.hh"
#include "../zerror_type/error_out.hh"

typedef struct{
	string name;
	SymmetricOperation data;
	ePointGroup gen_gp; // The group generated by the symmetric operation.
} symmetric_operation;

static const Int4 ISIZE_TABLE = 64;
static const symmetric_operation symmetric_operation_data[ISIZE_TABLE]={
		{ "Id",{{0,1,2,3},{1,1,1,1}},C1 },
		{ "Inv",{{0,1,2,3},{-1,-1,-1,-1}},Ci },
		{ "C2X",{{0,1,2,3},{1,-1,-1,1}},C2_X },
		{ "C2Y",{{0,1,2,3},{-1,1,-1,1}},C2_Y },
		{ "C2Z",{{0,1,2,3},{-1,-1,1,1}},C2_Z },
//		{ "C2",{{0,1,2,3},{-1,-1,1,-1}},C2_hex },
		{ "Sigma_X",{{0,1,2,3},{-1,1,1,1}},Cs_X },
		{ "Sigma_Y",{{0,1,2,3},{1,-1,1,1}},Cs_Y },
		{ "Sigma_Z",{{0,1,2,3},{1,1,-1,1}},Cs_Z },
//		{ "Sigma",{{0,1,2,3},{1,1,-1,1}},Cs_hex },
		{ "C4X+",{{0,2,1,3},{1,-1,1,1}},C4_X },
		{ "C4Y+",{{2,1,0,3},{1,1,-1,1}},C4_Y },
		{ "C4Z+",{{1,0,2,3},{-1,1,1,1}},C4_Z },
		{ "C4X-",{{0,2,1,3},{1,1,-1,1}},C4_X },
		{ "C4Y-",{{2,1,0,3},{-1,1,1,1}},C4_Y },
		{ "C4Z-",{{1,0,2,3},{1,-1,1,1}},C4_Z },
		{ "S4X-",{{0,2,1,3},{-1,1,-1,-1}},S4_X },
		{ "S4Y-",{{2,1,0,3},{-1,-1,1,-1}},S4_Y },
		{ "S4Z-",{{1,0,2,3},{1,-1,-1,-1}},S4_Z },
		{ "S4X+",{{0,2,1,3},{-1,-1,1,-1}},S4_X },
		{ "S4Y+",{{2,1,0,3},{1,-1,-1,-1}},S4_Y },
		{ "S4Z+",{{1,0,2,3},{-1,1,-1,-1}},S4_Z },
		{ "C31+",{{2,0,1,3},{1,1,1,1}},C31_rho },
		{ "C32+",{{2,0,1,3},{-1,1,-1,1}},C32_rho },
		{ "C33+",{{2,0,1,3},{-1,-1,1,1}},C33_rho },
		{ "C34+",{{2,0,1,3},{1,-1,-1,1}},C34_rho },
		{ "C31-",{{1,2,0,3},{1,1,1,1}},C31_rho },
		{ "C32-",{{1,2,0,3},{1,-1,-1,1}},C32_rho },
		{ "C33-",{{1,2,0,3},{-1,1,-1,1}},C33_rho },
		{ "C34-",{{1,2,0,3},{-1,-1,1,1}},C34_rho },
		{ "S61-",{{2,0,1,3},{-1,-1,-1,-1}},C31i_rho },
		{ "S62-",{{2,0,1,3},{1,-1,1,-1}},C32i_rho },
		{ "S63-",{{2,0,1,3},{1,1,-1,-1}},C33i_rho },
		{ "S64-",{{2,0,1,3},{-1,1,1,-1}},C34i_rho },
		{ "S61+",{{1,2,0,3},{-1,-1,-1,-1}},C31i_rho },
		{ "S62+",{{1,2,0,3},{-1,1,1,-1}},C32i_rho },
		{ "S63+",{{1,2,0,3},{1,-1,1,-1}},C33i_rho },
		{ "S64+",{{1,2,0,3},{1,1,-1,-1}},C34i_rho },
		{ "C6+",{{1,3,2,0},{-1,-1,1,-1}},C6 },
		{ "C3+",{{3,0,2,1},{1,1,1,1}},C3_hex },
		{ "C3-",{{1,3,2,0},{1,1,1,1}},C3_hex },
		{ "C6-",{{3,0,2,1},{-1,-1,1,-1}},C6 },
		{ "S3-",{{1,3,2,0},{1,1,-1,1}},C6h },
		{ "S6-",{{3,0,2,1},{-1,-1,-1,-1}},C3i_hex },
		{ "S6+",{{1,3,2,0},{-1,-1,-1,-1}},C3i_hex },
		{ "S3+",{{3,0,2,1},{1,1,-1,1}},C6h },
		{ "C2F",{{0,2,1,3},{-1,-1,-1,1}},C2D_X0_rho },
		{ "C2D",{{0,2,1,3},{-1,1,1,1}},C2D_X1_rho },
		{ "C2E",{{2,1,0,3},{-1,-1,-1,1}},C2D_Y0_rho },
		{ "C2C",{{2,1,0,3},{1,-1,1,1}},C2D_Y1_rho },
		{ "C2B(C'23)",{{1,0,2,3},{-1,-1,-1,1}},C2D_Z0 },
		{ "C2A(C\"23)",{{1,0,2,3},{1,1,-1,1}},C2D_Z1 },
		{ "SigmaF",{{0,2,1,3},{1,1,1,-1}},CsD_X0_rho },
		{ "SigmaD",{{0,2,1,3},{1,-1,-1,-1}},CsD_X1_rho },
		{ "SigmaE",{{2,1,0,3},{1,1,1,-1}},CsD_Y0_rho },
		{ "SigmaC",{{2,1,0,3},{-1,1,-1,-1}},CsD_Y1_rho },
		{ "SigmaB(SigmaD3)",{{1,0,2,3},{1,1,1,-1}},CsD_Z0 },
		{ "SigmaA(SigmaV3)",{{1,0,2,3},{-1,-1,1,-1}},CsD_Z1 },
		{ "C'21",{{0,3,2,1},{-1,-1,-1,-1}},C2D_X0_hex },
		{ "C\"21",{{0,3,2,1},{1,1,-1,1}},C2D_X1_hex },
		{ "C'22",{{3,1,2,0},{-1,-1,-1,-1}},C2D_Y0_hex },
		{ "C\"22",{{3,1,2,0},{1,1,-1,1}},C2D_Y1_hex },
//		{ "C'23",{{1,0,2,3},{-1,-1,-1,-1}},C2D_Z0_hex },
		{ "SigmaD1",{{0,3,2,1},{1,1,1,1}},CsD_X0_hex },
		{ "SigmaV1",{{0,3,2,1},{-1,-1,1,-1}},CsD_X1_hex },
		{ "SigmaD2",{{3,1,2,0},{1,1,1,1}},CsD_Y0_hex },
		{ "SigmaV2",{{3,1,2,0},{-1,-1,1,-1}},CsD_Y1_hex }
//		{ "SigmaD3",{{1,0,2,3},{1,1,1,1}},CsD_Z0_hex }
	};

const SymmetricOperation& change_enum_to_data(const eSymmetricOperation& num)
{
	if( num == D2dummy )
	{
        throw nerror_arg("D2dummy", __FILE__, __LINE__, __FUNCTION__);
	}
	return symmetric_operation_data[Int4(num)].data;
}

bool change_data_to_enum(const SymmetricOperation& symop, eSymmetricOperation& num)
{
	Int4 index;
	for(index=0; index<ISIZE_TABLE; index++)
		if( symmetric_operation_data[index].data == symop ){
			num = eSymmetricOperation(index);
			return true;
		}
	return false;
}

const string& Name(const eSymmetricOperation& num)
{
	if( num == D2dummy )
	{
        throw nerror_arg("D2dummy", __FILE__,  __LINE__, __FUNCTION__);
	}
	return symmetric_operation_data[Int4(num)].name;
}

const ePointGroup& generateGroup(const eSymmetricOperation& num)
{
	if( num == D2dummy )
	{
        throw nerror_arg("D2dummy", __FILE__, __LINE__, __FUNCTION__);
	}
	return symmetric_operation_data[Int4(num)].gen_gp;
}

void putAllElement(const ePointGroup& epg, vector<SymmetricOperation>& sym_opt)
{
	if( epg == C1 )
	{
		sym_opt.clear();
		sym_opt.push_back(change_enum_to_data(Id));
		return;
	}

	const eGroupToMaxSubgp& epg_to_sub_epg = enumMaxNormalSubgroup(epg);
	putAllElement(enumLowerGroup(epg_to_sub_epg), sym_opt);
	
	const Int4 isize = sym_opt.size();
	const Int4 isize2 = isize * 2;

	Int4 index;
	eSymmetricOperation esym;
	CosetRepresentativeMaxSubgp(epg_to_sub_epg, index,  esym);
	
	sym_opt.resize( isize * index );

	if( esym == D2dummy )
	{
		const SymmetricOperation sym_x = change_enum_to_data(C2X);
		for(Int4 k=0, k2=isize; k<isize; k++, k2++) sym_opt[k2] = sym_opt[k] * sym_x;  

		const SymmetricOperation sym_y = change_enum_to_data(C2Y);
		for(Int4 k=0, k2=isize2; k<isize2; k++, k2++) sym_opt[k2] = sym_opt[k] * sym_y;  
	}
	else
	{
		const SymmetricOperation sym = change_enum_to_data(esym);
		for(Int4 k=0, k2=isize; k<isize; k++, k2++) sym_opt[k2] = sym_opt[k] * sym;  
		if( index > 2 )
		{
			for(Int4 k=isize, k2=isize2; k<isize2; k++, k2++) sym_opt[k2] = sym_opt[k] * sym;  
		}
	}
}

void putSumAllElement(const ePointGroup& epg, NRMat<Int4>& ans)
{
	vector<SymmetricOperation> sym_opt;
	putAllElement(epg, sym_opt);
	
	ans = NRMat<Int4>(3, 3, 0);
	NRMat<Int4> mat(3, 3);
	for(vector<SymmetricOperation>::const_iterator it=sym_opt.begin(); it!=sym_opt.end(); it++)
	{
		putMatrixForm(*it, mat);
		for(Int4 i=0; i<3; i++)
			for(Int4 j=0; j<3; j++) ans[i][j] += mat[i][j];
	}
}

const SymmetricOperation& E()
{
	return change_enum_to_data(Id);
}

const SymmetricOperation& I()
{
	return change_enum_to_data(Inv);
}
