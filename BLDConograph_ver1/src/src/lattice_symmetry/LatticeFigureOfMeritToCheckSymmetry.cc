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
#include <limits>
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "check_equiv.hh"
#include "LatticeFigureOfMeritToCheckSymmetry.hh"

LatticeFigureOfMeritToCheckSymmetry::LatticeFigureOfMeritToCheckSymmetry()
	: m_S_red( SymMat43_Double( SymMat<Double>(3), NRMat<Int4>(4,3) ) ),
	  m_transform_to_lattice_equiv(3,3)
{
}


LatticeFigureOfMeritToCheckSymmetry::LatticeFigureOfMeritToCheckSymmetry(const BravaisType& brat,
		const SymMat43_Double& S)
	: m_S_red( SymMat43_Double( SymMat<Double>(3), NRMat<Int4>(4,3) ) ),
	  m_transform_to_lattice_equiv(3,3)
{
	setLatticeConstants43(brat, S);
	m_distance = numeric_limits<Double>::max();
}


void LatticeFigureOfMeritToCheckSymmetry::setDistance(const SymMat<Double>& S_super0)
{
	const NRMat<Int4> transMat = mprod( m_transform_to_lattice_equiv,
											put_transform_matrix_row4to3( m_S_red.second ) );
	const SymMat<Double> S_super_equiv = transform_sym_matrix(transMat, m_S_red.first);

	Double ans = 0.0;
	for(Int4 i=0; i<3; i++)
	{
		const Double diff = S_super0(i,i) - S_super_equiv(i,i);
		ans += diff * diff;

		for(Int4 j=0; j<i; j++)
		{
			const Double diff = S_super0(i,j) - S_super_equiv(i,j);
			ans += diff * diff * 2.0;
		}
	}
	m_distance = sqrt(ans);
}


#ifdef DEBUG
static bool checkInitialLatticeParameters(
		const BravaisType& brat,
		const SymMat43_Double& S_red)
{
	const SymMat<Double> dbl_S_red( S_red.first );

	if( brat.enumPointGroup() == Ci && brat.enumCentringType() == Prim )
	{
		assert( dbl_S_red(2,2)*0.9999 < dbl_S_red(1,1) && dbl_S_red(1,1)*0.9999 < dbl_S_red(0,0)
				&& fabs( dbl_S_red(0,1) ) * 1.9999 < dbl_S_red(1,1)
				&& fabs( dbl_S_red(0,2) ) * 1.9999 < dbl_S_red(2,2)
				&& fabs( dbl_S_red(1,2) ) * 1.9999 < dbl_S_red(2,2) );
	}
	else if( brat.enumPointGroup() == C2h_Y && brat.enumCentringType() == Prim )
	{
		assert( 0.0 <= dbl_S_red(0,2) && dbl_S_red(2,2)*0.9999 < dbl_S_red(0,0)
				&& fabs( dbl_S_red(0,2) ) * 1.9999 < dbl_S_red(2,2) && fabs( dbl_S_red(0,2) ) * 1.9999 < dbl_S_red(0,0) );
	}
	else if( brat.enumPointGroup() == C2h_Z && brat.enumCentringType() == Prim )
	{
		assert( 0.0 <= dbl_S_red(0,1) && dbl_S_red(1,1)*0.9999 < dbl_S_red(0,0)
				&& fabs( dbl_S_red(0,1) ) * 1.9999 < dbl_S_red(0,0) && fabs( dbl_S_red(0,1) ) * 1.9999 < dbl_S_red(1,1) );
	}
	else if( brat.enumPointGroup() == C2h_X && brat.enumCentringType() == Prim )
	{
		assert( 0.0 <= dbl_S_red(1,2) && dbl_S_red(2,2)*0.9999 < dbl_S_red(1,1)
				&& fabs( dbl_S_red(1,2) ) * 1.9999 < dbl_S_red(1,1) && fabs( dbl_S_red(1,2) ) * 1.9999 < dbl_S_red(2,2) );
	}
	else if( brat.enumPointGroup() == C2h_Y && brat.enumCentringType() == BaseZ )
	{
		assert( 0.0 <= dbl_S_red(0,2) && fabs( dbl_S_red(0,2) ) * 1.9999 < dbl_S_red(2,2) && fabs( dbl_S_red(0,2) ) * 0.9999 < dbl_S_red(0,0) );
	}
	else if( brat.enumPointGroup() == C2h_Z && brat.enumCentringType() == BaseX )
	{
		assert( 0.0 <= dbl_S_red(0,1) && fabs( dbl_S_red(0,1) ) * 1.9999 < dbl_S_red(0,0) && fabs( dbl_S_red(0,1) ) * 0.9999 < dbl_S_red(1,1) );
	}
	else if( brat.enumPointGroup() == C2h_X && brat.enumCentringType() == BaseY )
	{
		assert( 0.0 <= dbl_S_red(1,2) && fabs( dbl_S_red(1,2) ) * 1.9999 < dbl_S_red(1,1) && fabs( dbl_S_red(1,2) ) * 0.9999  < dbl_S_red(2,2) );
	}
	else if( brat.enumPointGroup() == D2h
			&& brat.enumCentringType() != BaseX
			&& brat.enumCentringType() != BaseY
			&& brat.enumCentringType() != BaseZ )
	{
		assert( dbl_S_red(2,2)*0.9999 < dbl_S_red(1,1) && dbl_S_red(1,1)*0.9999 < dbl_S_red(0,0) );
	}

	const SymMat<Double> S_super = transform_sym_matrix(S_red.second, S_red.first);
	assert( S_super(0,1) <= 0.0
			&& S_super(0,2) <= 0.0
			&& S_super(0,3) <= 0.0
			&& S_super(1,2) <= 0.0
			&& S_super(1,3) <= 0.0
			&& S_super(2,3) <= 0.0 );

	return true;
}
#endif


void LatticeFigureOfMeritToCheckSymmetry::setLatticeConstants43(const BravaisType& brat, const SymMat43_Double& S)
{
	m_S_red = S;
	assert( checkInitialLatticeParameters(brat, m_S_red) );

	m_latfom.setLatticeConstants43(brat, S);
	m_distance = numeric_limits<Double>::max();
	m_transform_to_lattice_equiv = 0;
}




bool LatticeFigureOfMeritToCheckSymmetry::checkIfLatticeIsMonoclinic(const ePointGroup& epg_new,
		const Double& resol,
		map< SymMat<Double>, NRMat<Int4> >& ans) const
{
	ans.clear();

	SymMat<Double> ans0 = m_S_red.first;
	cal_average_crystal_system(C2h_X, ans0);

	SymMat<Double> S_red(3);
	NRMat<Int4> trans_mat2;
	if( check_equiv_m(ans0, m_S_red.first, resol ) )
	{
		if( epg_new == C2h_X )
		{
			S_red = ans0;
			trans_mat2 = m_S_red.second;
			putBuergerReducedMonoclinicP(1, 2, S_red, trans_mat2);
		}
		else if( epg_new == C2h_Y )
		{
			S_red = transform_sym_matrix(put_matrix_YXZ(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_YXZ());
			putBuergerReducedMonoclinicP(0, 2, S_red, trans_mat2);
		}
		else // if( epg_new == C2h_Z )
		{
			S_red = transform_sym_matrix(put_matrix_YZX(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_ZXY());
			putBuergerReducedMonoclinicP(0, 1, S_red, trans_mat2);
		}
		ans.insert( SymMat43_Double( S_red, trans_mat2) );
	}

	ans0 = m_S_red.first;
	cal_average_crystal_system(C2h_Y, ans0);
	if( check_equiv_m(ans0, m_S_red.first, resol ) )
	{
		if( epg_new == C2h_X )
		{
			S_red = transform_sym_matrix(put_matrix_YXZ(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_YXZ());
			putBuergerReducedMonoclinicP(1, 2, S_red, trans_mat2);
		}
		else if( epg_new == C2h_Y )
		{
			S_red = ans0;
			trans_mat2 = m_S_red.second;
			putBuergerReducedMonoclinicP(0, 2, S_red, trans_mat2);
		}
		else // if( epg_new == C2h_Z )
		{
			S_red = transform_sym_matrix(put_matrix_XZY(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_XZY());
			putBuergerReducedMonoclinicP(0, 1, S_red, trans_mat2);
		}
		ans.insert( SymMat43_Double( S_red, trans_mat2) );
	}

	ans0 = m_S_red.first;
	cal_average_crystal_system(C2h_Z, ans0);
	if( check_equiv_m(ans0, m_S_red.first, resol ) )
	{
		if( epg_new == C2h_X )
		{
			S_red = transform_sym_matrix(put_matrix_ZXY(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_YZX());
			putBuergerReducedMonoclinicP(1, 2, S_red, trans_mat2);
		}
		else if( epg_new == C2h_Y )
		{
			S_red = transform_sym_matrix(put_matrix_XZY(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_XZY());
			putBuergerReducedMonoclinicP(0, 2, S_red, trans_mat2);
		}
		else // if( epg_new == C2h_Z )
		{
			S_red = ans0;
			trans_mat2 = m_S_red.second;
			putBuergerReducedMonoclinicP(0, 1, S_red, trans_mat2);
		}
		ans.insert( SymMat43_Double( S_red, trans_mat2) );
	}

	return !( ans.empty() );
}


bool LatticeFigureOfMeritToCheckSymmetry::checkIfLatticeIsOrthorhombic(const Double& resol,
							map< SymMat<Double>, NRMat<Int4> >& ans) const
{
	ans.clear();

	const BravaisType& brat = m_latfom.putBravaisType();

	SymMat<Double> ans0 = m_S_red.first;
	cal_average_crystal_system(D2h, ans0);
	if( check_equiv_m(ans0, m_S_red.first, resol ) )
	{
		if( brat.enumCentringType() == BaseX )
		{
			if( ans0(1,1) < ans0(2,2) )
			{
				ans.insert( SymMat43_Double( transform_sym_matrix(put_matrix_ZYX(), ans0), mprod( m_S_red.second, put_matrix_ZYX() ) ) );
			}
			else
			{
				ans.insert( SymMat43_Double( transform_sym_matrix(put_matrix_YZX(), ans0), mprod( m_S_red.second, put_matrix_ZXY() ) ) );
			}
		}
		else if( brat.enumCentringType() == BaseY )
		{
			if( ans0(0,0) < ans0(2,2) )
			{
				ans.insert( SymMat43_Double( transform_sym_matrix(put_matrix_ZXY(), ans0), mprod( m_S_red.second, put_matrix_YZX() ) ) );
			}
			else
			{
				ans.insert( SymMat43_Double( transform_sym_matrix(put_matrix_XZY(), ans0), mprod( m_S_red.second, put_matrix_XZY() ) ) );
			}
		}
		else if( brat.enumCentringType() == BaseZ )
		{
			if( ans0(0,0) < ans0(1,1) )
			{
				ans.insert( SymMat43_Double( transform_sym_matrix(put_matrix_YXZ(), ans0), mprod( m_S_red.second, put_matrix_YXZ() ) ) );
			}
			else
			{
				ans.insert( SymMat43_Double( transform_sym_matrix(put_matrix_XYZ(), ans0), m_S_red.second ) );
			}
		}
		else
		{
			NRMat<Int4> trans_mat = m_S_red.second;
			putBuergerReducedOrthorhombic(brat.enumCentringType(), ans0, trans_mat);
			ans.insert( SymMat43_Double(ans0, trans_mat ) );
		}
		return true;
	}
	return false;
}


bool LatticeFigureOfMeritToCheckSymmetry::checkIfLatticeIsTetragonal(const Double& resol,
		map< SymMat<Double>, NRMat<Int4> >& ans) const
{
	ans.clear();

	SymMat<Double> ans0 = m_S_red.first;
	cal_average_crystal_system(D4h_X, ans0);
	if( check_equiv_m(ans0, m_S_red.first, resol ) )
	{
		ans.insert( SymMat43_Double(
				transform_sym_matrix(put_matrix_YZX(), ans0), mprod( m_S_red.second, put_matrix_ZXY() ) ) );
	}

	ans0 = m_S_red.first;
	cal_average_crystal_system(D4h_Y, ans0);
	if( check_equiv_m(ans0, m_S_red.first, resol ) )
	{
		ans.insert( SymMat43_Double(
				transform_sym_matrix(put_matrix_XZY(), ans0), mprod( m_S_red.second, put_matrix_XZY() ) ) );
	}

	ans0 = m_S_red.first;
	cal_average_crystal_system(D4h_Z, ans0);
	if( check_equiv_m(ans0, m_S_red.first, resol ) )
	{
		ans.insert( SymMat43_Double(ans0, m_S_red.second ) );
	}

	return !( ans.empty() );
}




bool LatticeFigureOfMeritToCheckSymmetry::checkIfLatticeIsHexagonal(const ePointGroup& epg_new, const Double& resol,
		map< SymMat<Double>, NRMat<Int4> >& ans) const
{
	ans.clear();
	const BravaisType& brat = m_latfom.putBravaisType();

	SymMat43_Double ans2(SymMat<Double>(3), NRMat<Int4>(3,3));

	if( brat.enumPointGroup() == C2h_X )
	{
		ans2.first = transform_sym_matrix(put_matrix_YZX(), m_S_red.first);
		ans2.second = mprod( m_S_red.second, put_matrix_ZXY() );
	}
	else if( brat.enumPointGroup() == C2h_Y )
	{
		ans2.first = transform_sym_matrix(put_matrix_ZXY(), m_S_red.first);
		ans2.second = mprod( m_S_red.second, put_matrix_YZX() );
	}
	else // if( brat.enumPointGroup() == C2h_Z )
	{
		ans2.first = transform_sym_matrix(put_matrix_XYZ(), m_S_red.first);
		ans2.second = m_S_red.second;
	}

	if( ans2.first(0,1) < 0.0 )
	{
		ans2.first(0,1) *= -1;
		ans2.second[0][0] *= -1;
		ans2.second[1][0] *= -1;
		ans2.second[2][0] *= -1;
	}

	SymMat<Double> ans0 = ans2.first;
	cal_average_crystal_system(epg_new, ans2.first);
	if( check_equiv_m(ans2.first, ans0, resol ) )
	{
		ans.insert( ans2 );
		return true;
	}
	else return false;
}


bool LatticeFigureOfMeritToCheckSymmetry::checkLatticeSymmetry(const ePointGroup& epg_new, const Double& resol,
		map< SymMat<Double>, NRMat<Int4> >& ans) const
{
	ans.clear();
	const BravaisType& brat = m_latfom.putBravaisType();
	if( epg_new == Ci || epg_new == brat.enumPointGroup() )
	{
		ans.insert( m_S_red );
		return true;
	}

	if( epg_new == C2h_X || epg_new == C2h_Y ||  epg_new == C2h_Z )
	{
		assert( brat.enumPointGroup() == Ci );
		assert( brat.enumCentringType() == Prim );

		return checkIfLatticeIsMonoclinic(epg_new, resol, ans);
	}
	else if( epg_new == D4h_Z )
	{
		assert( brat.enumPointGroup() == D2h );
		assert( brat.enumCentringType() == Prim
				|| brat.enumCentringType() == Inner );

		return checkIfLatticeIsTetragonal(resol, ans);
	}
	else if( epg_new == D2h )
	{
		assert( brat.enumPointGroup() != Ci || brat.enumCentringType() == Prim );
		assert( brat.enumPointGroup() != C2h_Z || brat.enumCentringType() == BaseX );
		assert( brat.enumPointGroup() != C2h_X || brat.enumCentringType() == BaseY );
		assert( brat.enumPointGroup() != C2h_Y || brat.enumCentringType() == BaseZ );
		assert( brat.enumCentringType() != Rhom_hex );

		return checkIfLatticeIsOrthorhombic(resol, ans);
	}
	else if( epg_new == D6h )
	{
		assert( brat.enumCentringType() == Prim );
		assert( brat.enumPointGroup() == C2h_X
				|| brat.enumPointGroup() == C2h_Y
				|| brat.enumPointGroup() == C2h_Z );
		return checkIfLatticeIsHexagonal(epg_new, resol, ans);
	}
	else
	{
		assert( epg_new == Oh );
		assert( brat.enumCentringType() == Prim
				|| brat.enumCentringType() == Inner
				|| brat.enumCentringType() == Face );

		SymMat43_Double ans2 = m_S_red;
		cal_average_crystal_system(epg_new, ans2.first);
		if( check_equiv_m(ans2.first, m_S_red.first, resol ) )
		{
			ans.insert( ans2 );
			return true;
		}
	}
	return !(ans.empty());
}


void LatticeFigureOfMeritToCheckSymmetry::putLatticesOfHigherSymmetry(
		const ePointGroup& epg, const Double& resol,
		vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result) const
{
	lattice_result.clear();
	map< SymMat<Double>, NRMat<Int4> > S_red_tray;
	if( !this->checkLatticeSymmetry(epg, resol, S_red_tray) ) return;

	const BravaisType& ebrat_original = this->putLatticeFigureOfMerit().putBravaisType();
	const eCentringType eblat = (ebrat_original.enumBravaisType()==Monoclinic_B?
									(epg==D31d_rho?Prim:(epg==D3d_1_hex?Rhom_hex:BaseZ)):ebrat_original.enumCentringType());

	const NRMat<Int4> matrix_min_to_sell = this->putInitialForm().second;

	SymMat<Double> S_super(4);
	NRMat<Int4> trans_mat(4,3);

	for(map< SymMat<Double>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
	{
		// S_super = it->second * it->first * Transpose(it->second) is close to
		// Delone-reduced form of the original lattice.
		S_super = transform_sym_matrix(it->second, it->first );

		trans_mat = identity_matrix<Int4>(4);

		// S_super = trans_mat * it->second * it->first * Transpose(trans_mat * it->second).
		put_Selling_reduced_dim_3(S_super, trans_mat);
		moveSmallerDiagonalLeftUpper(S_super, trans_mat);

		lattice_result.push_back( LatticeFigureOfMeritToCheckSymmetry( BravaisType( pair<eCentringType, ePointGroup>(eblat, epg) ),
										SymMat43_Double(it->first, mprod(trans_mat, it->second) ) ) );

		lattice_result.rbegin()->setTransformToOriginalLattice( mprod(this->putTransformToOriginalLattice(), Inverse3( put_transform_matrix_size4to3(trans_mat) ) ) );
	}
}


void LatticeFigureOfMeritToCheckSymmetry::printLatticeInformation(
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis,
					const Int4& label_start0,
					ostream* os) const
{
	m_latfom.printLatticeInformation(abc_axis, rh_axis, label_start0, os);


  	Int4 label_start = label_start0;
	os->width(label_start);
	*os << "" << "<DistanceFromInputUnitCell>";
  	os->width(14);
   	*os << sqrt(m_distance);
  	*os << " </DistanceFromInputUnitCell>\n";
}
