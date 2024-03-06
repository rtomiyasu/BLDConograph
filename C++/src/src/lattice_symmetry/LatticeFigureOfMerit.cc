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
#include "../utility_data_structure/FracMat.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "ReducedLatticeToCheckBravais.hh"
#include "LatticeFigureOfMerit.hh"

const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_face = put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToFace() ) );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_body = put_transform_matrix_row3to4( BravaisType::putTransformMatrixFromBodyToPrimitive() );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_rhomhex = put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToRhomHex() ) );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_base[3] =
		{
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseA_Axis) ) ),
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseB_Axis) ) ),
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseC_Axis) ) )
		};
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_prim = put_transform_matrix_row3to4();

LatticeFigureOfMerit::LatticeFigureOfMerit()
	: m_S_red(3)
{
}


LatticeFigureOfMerit::LatticeFigureOfMerit(const BravaisType& brat,
		const SymMat43_Double& S)
	: m_S_red(3)
{
	this->setLatticeConstants43(brat, S);
}

#ifdef DEBUG
static bool checkInitialLatticeParameters(
		const BravaisType& brat,
		const SymMat<Double>& S_red)
{
	const SymMat<Double> inv_S_red( Inverse3(S_red) );

	if( brat.enumPointGroup() == C2h_Y && brat.enumCentringType() == Prim )
	{
		assert( inv_S_red(0,2) <= 0.0 &&
				inv_S_red(0,0) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(0,0) );
	}
	else if( brat.enumPointGroup() == C2h_Z && brat.enumCentringType() == Prim )
	{
		assert( inv_S_red(0,1) <= 0.0
				&& inv_S_red(0,0) * 0.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(0,0)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumPointGroup() == C2h_X && brat.enumCentringType() == Prim )
	{
		assert( inv_S_red(1,2) <= 0.0
				&& inv_S_red(1,1) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(2,2) );
	}
	else if( brat.enumPointGroup() == C2h_Y && brat.enumCentringType() == BaseZ )
	{
		assert( inv_S_red(0,2) <= 0.0
				&& fabs( inv_S_red(0,2) ) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(0,0) );
	}
	else if( brat.enumPointGroup() == C2h_Z && brat.enumCentringType() == BaseX )
	{
		assert( inv_S_red(0,1) <= 0.0
				&& fabs( inv_S_red(0,1) ) * 0.9999 < inv_S_red(0,0)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumPointGroup() == C2h_X && brat.enumCentringType() == BaseY )
	{
		assert( inv_S_red(1,2) <= 0.0
				&& fabs( inv_S_red(1,2) ) * 0.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(2,2) );
	}
	else if( brat.enumBravaisType() == Orthorhombic_C )
	{
		assert( brat.enumCentringType() == BaseZ );
		assert( inv_S_red(0,0) * 0.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumPointGroup() == D2h && brat.enumCentringType() == Prim )
	{
		assert( inv_S_red(0,0) * 0.9999 < inv_S_red(1,1)
				&& inv_S_red(1,1) * 0.9999 < inv_S_red(2,2) );
	}
	return true;
}
#endif

static void putTransformMatrixToBuergerReduced(
		const SymMat<Double>& S, NRMat<Int4>& trans_mat)
{
	assert( S.size() == 3 );

	SymMat<Double> S_super_obtuse(4);
	put_Selling_reduced_dim_3(S, S_super_obtuse, trans_mat);
	moveSmallerDiagonalLeftUpper(S_super_obtuse, trans_mat);

	// S_red = trans_mat * S_super_obtuse * transpose(trans_mat).
	SymMat<Double> S_red(3);
	NRMat<Int4> trans_mat2;
	putBuergerReducedMatrix(S_super_obtuse, S_red, trans_mat2);
	trans_mat = mprod( trans_mat2, put_transform_matrix_row4to3(trans_mat) );
}


void LatticeFigureOfMerit::setInverseOfBuergerReducedForm(NRMat<Int4>& trans_mat, const SymMat43_Double& S_optimized)
{
	if( m_brat.enumBravaisType() == Triclinic )
	{
		// trans_mat * Inverse(S_optimized.first) * transpose(trans_mat) is Buerger-reduced
		// <=> Inverse of transpose(Inverse(trans_mat)) * S_optimized.first * Inverse(trans_mat) is Buerger-reduced.
		putTransformMatrixToBuergerReduced(Inverse3(S_optimized.first), trans_mat);
		transpose_square_matrix(trans_mat);
		m_S_red = transform_sym_matrix(Inverse3(trans_mat), S_optimized.first);
	}
	else
	{
		m_S_red = S_optimized.first;
		trans_mat = identity_matrix<Int4>(3);
		if( m_brat.enumBravaisType() == Monoclinic_P )
		{
			if( m_brat.enumPointGroup() == C2h_X )
			{
				putBuergerReducedMonoclinicP(1, 2, m_S_red, trans_mat);
			}
			else if( m_brat.enumPointGroup() == C2h_Y )
			{
				putBuergerReducedMonoclinicP(0, 2, m_S_red, trans_mat);
			}
			else //if( m_brat.enumPointGroup() == C2h_Z )
			{
				putBuergerReducedMonoclinicP(0, 1, m_S_red, trans_mat);
			}
		}
		else if( m_brat.enumBravaisType() == Monoclinic_B )
		{
			m_S_red = S_optimized.first;
			putBuergerReducedMonoclinicB(m_brat, m_S_red, trans_mat);
		}
		else if( m_brat.enumPointGroup() == D2h )
		{
			m_S_red = S_optimized.first;
			putBuergerReducedOrthorhombic(m_brat.enumCentringType(), m_S_red, trans_mat);
		}
	}

	assert( checkInitialLatticeParameters(m_brat, m_S_red) );
}


void LatticeFigureOfMerit::setLatticeConstants43(const BravaisType& brat, const SymMat43_Double& S)
{
	m_brat = brat;

	NRMat<Int4> trans_mat;
	setInverseOfBuergerReducedForm(trans_mat, S);	// Set m_S_red from S.
}

ZErrorMessage LatticeFigureOfMerit::setLatticeConstants(const BravaisType& brat, const SymMat<Double>& Sval)
{
	assert( Sval.size()==3 );

	SymMat43_Double S_red_optimized = SymMat43_Double(Sval, NRMat<Int4>(4,3));
	cal_average_crystal_system(brat.enumPointGroup(), S_red_optimized.first);
	if( brat.enumCentringType() == Face )
	{
		S_red_optimized.second = m_tmat_prim_to_face;
	}
	else if( brat.enumCentringType() == Inner )
	{
		S_red_optimized.second = m_tmat_prim_to_body;
	}
	else if( brat.enumCentringType() == BaseX
			|| brat.enumCentringType() == BaseY
			|| brat.enumCentringType() == BaseZ )
	{
		S_red_optimized.second = m_tmat_prim_to_base[ (size_t)brat.enumBASEaxis() ];
	}
	else if( brat.enumCentringType() == Rhom_hex )
	{
		S_red_optimized.second = m_tmat_prim_to_rhomhex;
	}
	else // if( brat.enumCentringType() == Prim )
	{
		S_red_optimized.second = m_tmat_prim_to_prim;
	}

	// S_super_obtuse = trans_mat * S_red.first * Transpose(trans_mat).
	SymMat<Double> S_super_obtuse = transform_sym_matrix(S_red_optimized.second, S_red_optimized.first);
	if( !put_Selling_reduced_dim_3(S_super_obtuse, S_red_optimized.second) )
	{
		return ZErrorMessage(ZErrorArgument, "The argument matrix is not positive definite" __FILE__, __LINE__, __FUNCTION__);
	}
	moveSmallerDiagonalLeftUpper(S_super_obtuse, S_red_optimized.second);

	setLatticeConstants43(brat, S_red_optimized);

	return ZErrorMessage();
}


inline bool checkIfFirstEntryIsPositive(const VecDat3<Int4>& rhs)
{
	for(Int4 i=0; i<3; i++)
	{
		if( rhs[i] == 0 ) continue;
		if( rhs[i] > 0 ) return true;
		else return false;
	}
	return false;
}



void LatticeFigureOfMerit::printLatticeInformation(
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis,
					const Int4& label_start0,
					ostream* os) const
{
	Int4 label_start = label_start0;
	os->width(label_start);
  	*os << "" << "<CrystalSystem>";
	os->width(17);
	*os << put_bravais_type_name(this->enumBravaisType(), abc_axis);
  	*os << " </CrystalSystem>\n\n";

	os->width(label_start); *os << "";
  	*os << "<!-- a, b, c(angstrom), alpha, beta, gamma(deg.)-->\n";

	VecDat3<Double> length_axis, angle_axis;
	if( this->enumBravaisType() == Rhombohedral )
	{
		this->putReducedLatticeConstantsDegree(abc_axis, Rho_Axis, length_axis, angle_axis);

		os->width(label_start); *os << "";
		*os << "<UnitCellParameters axis=\"Rhombohedral\">";
	  	os->width(14);
	  	*os << length_axis[0];
	  	os->width(14);
	   	*os << length_axis[1];
	  	os->width(14);
	   	*os << length_axis[2];
	 	os->width(14);
	   	*os << angle_axis[0];
	  	os->width(14);
	   	*os << angle_axis[1];
	  	os->width(14);
	   	*os << angle_axis [2];
	  	*os << " </UnitCellParameters>\n";

		this->putReducedLatticeConstantsDegree(abc_axis, Hex_Axis, length_axis, angle_axis);

		os->width(label_start); *os << "";
	  	*os << "<UnitCellParameters axis=\"Hexagonal\">";
	  	os->width(14);
	  	*os << length_axis[0];
	  	os->width(14);
	   	*os << length_axis[1];
	  	os->width(14);
	   	*os << length_axis[2];
	 	os->width(14);
	   	*os << angle_axis[0];
	  	os->width(14);
	   	*os << angle_axis[1];
	  	os->width(14);
	   	*os << angle_axis[2];
	  	*os << " </UnitCellParameters>\n\n";
	}
	else
	{
		this->putReducedLatticeConstantsDegree(abc_axis, Rho_Axis, length_axis, angle_axis);

		os->width(label_start); *os << "";
		*os << "<UnitCellParameters>";
	  	os->width(14);
	  	*os << length_axis[0];
	  	os->width(14);
	   	*os << length_axis[1];
	  	os->width(14);
	   	*os << length_axis[2];
	 	os->width(14);
	   	*os << angle_axis[0];
	  	os->width(14);
	   	*os << angle_axis[1];
	  	os->width(14);
	   	*os << angle_axis[2];
	  	*os << " </UnitCellParameters>\n";
	}
}


void LatticeFigureOfMerit::putLatticeConstantsDegree(const BravaisType& brat, const SymMat<Double>& S0,
		const eABCaxis& axis1,
		const eRHaxis& axis2, VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis)
{
	SymMat<Double> S = S0;
	if( brat.enumBravaisType() == Rhombohedral && axis2 != brat.enumRHaxis() )
	{
		if( axis2 == Hex_Axis ) // Rho -> Hex.
		{
			static const FracMat matrix_rho2hex = FInverse3( transpose(BravaisType::putTransformMatrixFromPrimitiveToRhomHex() ) );
			S = transform_sym_matrix(matrix_rho2hex.mat, S)/(matrix_rho2hex.denom*matrix_rho2hex.denom);
		}
		else // if( axis2 == RhoAxis ) // Hex -> Rho.
		{
			static const NRMat<Int4> matrix_hex2rho = transpose( BravaisType::putTransformMatrixFromPrimitiveToRhomHex() );
			S = transform_sym_matrix(matrix_hex2rho, S);
		}
	}
	else if( brat.enumBravaisType() == Monoclinic_B )
	{
		const NRMat<Int4> this2output = put_transform_matrix_monoclinic_b(brat.enumABCaxis(), axis1);
		S = transform_sym_matrix(this2output, S);
	}

	calLatticeConstant( S, length_axis, angle_axis );
}
