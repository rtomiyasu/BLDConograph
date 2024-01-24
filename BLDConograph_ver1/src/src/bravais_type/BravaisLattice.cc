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
#include "../utility_data_structure/index_set.hh"
#include "../utility_func/transform_sym_matrix.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "../lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.hh"
#include "../lattice_symmetry/ReducedLatticeToCheckBravais.hh"
#include "../zerror_type/error_out.hh"
#include "../zlog/zlog.hh"
#include "../ControlParam.hh"
#include "BravaisLattice.hh"

const bool BravaisLattice::m_DoesPrudentSymSearch = false;

BravaisLattice::BravaisLattice()
{
	for(Int4 i=0; i<NUM_LS; i++)
	{
		OutputSymmetry[i] = false;
		JudgeSymmetry[i] = false;
	}
	
   	m_resol = 0.0;
}


BravaisLattice::~BravaisLattice()
{
}


// Set the member variables.
void BravaisLattice::setParam(const ControlParam& cont) 
{
	OutputSymmetry[(size_t)Triclinic] = cont.putOutputSymmetry(Triclinic);
	JudgeSymmetry[(size_t)Triclinic] = false;
	for(Int4 i=1; i<NUM_LS; i++)
	{
		OutputSymmetry[i] = cont.putOutputSymmetry(eBravaisType(i));
		JudgeSymmetry[i] = cont.putOutputSymmetry(eBravaisType(i));
	}

	if( JudgeSymmetry[(size_t)Cubic_P] )
	{
		JudgeSymmetry[(size_t)Tetragonal_P] = true;
	}
	if( JudgeSymmetry[(size_t)Hexagonal] )
	{
		JudgeSymmetry[(size_t)Monoclinic_P] = true;
	}
	if( JudgeSymmetry[(size_t)Tetragonal_P] )
	{
		JudgeSymmetry[(size_t)Orthorhombic_P] = true;
	}
	if( JudgeSymmetry[(size_t)Orthorhombic_P] )
	{
		JudgeSymmetry[(size_t)Monoclinic_P] = true;
	}
	
	if( JudgeSymmetry[(size_t)Orthorhombic_C] )
	{
		JudgeSymmetry[(size_t)Monoclinic_B] = true;
	}

	if( JudgeSymmetry[(size_t)Cubic_I] )
	{
		JudgeSymmetry[(size_t)Tetragonal_I] = true;
	}
	if( JudgeSymmetry[(size_t)Tetragonal_I] )
	{
		JudgeSymmetry[(size_t)Orthorhombic_I] = true;
	}

	if( JudgeSymmetry[(size_t)Cubic_F] )
	{
		JudgeSymmetry[(size_t)Orthorhombic_F] = true;
	}
	
	m_resol = cont.putResolution();
}




void BravaisLattice::putCentringTypes(const ReducedLatticeToCheckBravais& RLCB,
		const BravaisType& brat,
		vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result)
{
	lattice_result.clear();
	
	const map< SymMat<Double>, NRMat<Int4> >& S_red_tray = RLCB.checkCentringType(brat);
	if( S_red_tray.empty() ) return;

	// The lattice of RLCB has at least the symmetry given by eblat.
	SymMat<Double> S_super(4);
	NRMat<Int4> trans_mat(4,3);

	for(map< SymMat<Double>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
	{
		S_super = transform_sym_matrix(it->second, it->first);
		trans_mat = identity_matrix<Int4>(4);

		// S_super = trans_mat * it->second * it->first * Transpose(trans_mat * it->second) is Delone reduced.
		if( !put_Selling_reduced_dim_3(S_super, trans_mat) )
		{
			assert( false );
		}
		moveSmallerDiagonalLeftUpper(S_super, trans_mat);

		lattice_result.push_back( LatticeFigureOfMeritToCheckSymmetry( brat, SymMat43_Double(it->first, mprod(trans_mat, it->second) ) ) );

		lattice_result.rbegin()->setTransformToOriginalLattice( Inverse3( put_transform_matrix_size4to3(trans_mat) ) );
	}
}



static SymMat43_Double putTransformMatrixFromSellingReducedToBuergerReduced(
		const SymMat<Double>& S_super_obtuse)
{
	NRMat<Int4> trans_mat(4,3);

	// S_red = trans_mat * S_super_obtuse * transpose(trans_mat).
	SymMat<Double> S_red(3);
	putBuergerReducedMatrix(S_super_obtuse, false, S_red, trans_mat);

	return SymMat43_Double( S_red, put_transform_matrix_row3to4( Inverse3( trans_mat ) ) );
}



void BravaisLattice::putLatticeCandidatesForEachBravaisTypes(const SymMat<Double>& S_super,
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis,
					vector<LatticeFigureOfMeritToCheckSymmetry> lattice_result[NUM_LS]) const
{
	try{

	for(Int4 i=1; i<NUM_LS; i++)
	{
		lattice_result[i].clear();
	}
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tri = lattice_result[(size_t)Triclinic];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_mono_P = lattice_result[(size_t)Monoclinic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_mono_B = lattice_result[(size_t)Monoclinic_B];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_P = lattice_result[(size_t)Orthorhombic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_B = lattice_result[(size_t)Orthorhombic_C];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_I = lattice_result[(size_t)Orthorhombic_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_F = lattice_result[(size_t)Orthorhombic_F];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tetra_P = lattice_result[(size_t)Tetragonal_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tetra_I = lattice_result[(size_t)Tetragonal_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_rhom = lattice_result[(size_t)Rhombohedral];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_hex = lattice_result[(size_t)Hexagonal];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_P = lattice_result[(size_t)Cubic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_I = lattice_result[(size_t)Cubic_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_F = lattice_result[(size_t)Cubic_F];

	const SymMat43_Double S_red = putTransformMatrixFromSellingReducedToBuergerReduced(S_super);
	LatticeFigureOfMeritToCheckSymmetry latFOM(BravaisType( pair<eCentringType, ePointGroup>(Prim, Ci) ), S_red);

	lattice_result_tri.push_back( latFOM );

	// Calculate figures of merit as triclinic
	const ReducedLatticeToCheckBravais RLCB(abc_axis, rh_axis, m_DoesPrudentSymSearch, m_resol, S_red);

	vector<LatticeFigureOfMeritToCheckSymmetry> latFOM_tray;

	if( JudgeSymmetry[Monoclinic_B] )
	{
		putCentringTypes(RLCB, BravaisType( put_monoclinic_b_type(abc_axis) ), latFOM_tray);
		lattice_result_mono_B.insert(lattice_result_mono_B.end(), latFOM_tray.begin(), latFOM_tray.end());
	}
	if( JudgeSymmetry[Orthorhombic_I] )
	{
		putCentringTypes(RLCB, BravaisType( pair<eCentringType, ePointGroup>(Inner, D2h) ), latFOM_tray);
		lattice_result_ortho_I.insert(lattice_result_ortho_I.end(), latFOM_tray.begin(), latFOM_tray.end());
	}
	if( JudgeSymmetry[Orthorhombic_F] )
	{
		putCentringTypes(RLCB, BravaisType( pair<eCentringType, ePointGroup>(Face, D2h) ), latFOM_tray);
		lattice_result_ortho_F.insert(lattice_result_ortho_F.end(), latFOM_tray.begin(), latFOM_tray.end());
	}
	if( JudgeSymmetry[Rhombohedral] )
	{
		putCentringTypes(RLCB, BravaisType( put_rhombohedral_type(rh_axis) ), latFOM_tray);
		lattice_result_rhom.insert(lattice_result_rhom.end(), latFOM_tray.begin(), latFOM_tray.end());
	}
	if( JudgeSymmetry[Monoclinic_P] )
	{
		latFOM.putLatticesOfHigherSymmetry(put_monoclinic_p_type(abc_axis), m_resol, latFOM_tray);
		lattice_result_mono_P.insert(lattice_result_mono_P.end(), latFOM_tray.begin(), latFOM_tray.end());
	}
	if( JudgeSymmetry[Orthorhombic_P] )
	{
		latFOM.putLatticesOfHigherSymmetry(D2h, m_resol, latFOM_tray);
		lattice_result_ortho_P.insert(lattice_result_ortho_P.end(), latFOM_tray.begin(), latFOM_tray.end());
	}

	for(Int4 i=1; i<NUM_LS; i++)
	{
		if( !JudgeSymmetry[i] ) continue;

		const Int4 num_lattice = lattice_result[i].size();

		for(Int4 index=0; index<num_lattice; index++)
		{
			LatticeFigureOfMeritToCheckSymmetry& latFOM0 = lattice_result[i][index];

			if( eBravaisType(i) == Monoclinic_P )
			{
				if( JudgeSymmetry[Hexagonal] )
				{
					latFOM0.putLatticesOfHigherSymmetry(D6h, m_resol, latFOM_tray);
					lattice_result_hex.insert(lattice_result_hex.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
			else if( eBravaisType(i) == Monoclinic_B )
			{
				if( JudgeSymmetry[Orthorhombic_C] )
				{
					latFOM0.putLatticesOfHigherSymmetry(D2h, m_resol, latFOM_tray);
					lattice_result_ortho_B.insert(lattice_result_ortho_B.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
			else if( eBravaisType(i) == Orthorhombic_P )
			{
				if( JudgeSymmetry[Tetragonal_P] )
				{
					latFOM0.putLatticesOfHigherSymmetry(D4h_Z, m_resol, latFOM_tray);
					lattice_result_tetra_P.insert(lattice_result_tetra_P.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
				if( JudgeSymmetry[Cubic_P] )
				{
					latFOM0.putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray);
					lattice_result_cubic_P.insert(lattice_result_cubic_P.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
			else if( eBravaisType(i) == Orthorhombic_I )
			{
				if( JudgeSymmetry[Tetragonal_I] )
				{
					latFOM0.putLatticesOfHigherSymmetry(D4h_Z, m_resol, latFOM_tray);
					lattice_result_tetra_I.insert(lattice_result_tetra_I.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
				if( JudgeSymmetry[Cubic_I] )
				{
					latFOM0.putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray);
					lattice_result_cubic_I.insert(lattice_result_cubic_I.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
			else if( eBravaisType(i) == Orthorhombic_F )
			{
				if( JudgeSymmetry[Cubic_F] )
				{
					latFOM0.putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray);
					lattice_result_cubic_F.insert(lattice_result_cubic_F.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
		}

ZLOG_INFO( "(" + num2str( i+1 ) + ") The number of candidates for " + put_bravais_type_name(eBravaisType(i), abc_axis)
			+ " : " + num2str<Int4>( lattice_result[i].size() ) + "\n" );
	}
ZLOG_INFO( "\n" );
    }
    catch(bad_alloc& ball){
    	throw nerror(ball, __FILE__, __LINE__, __FUNCTION__);
    }
}


void BravaisLattice::execute(const SymMat<Double>& S_super,
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis,
					vector<LatticeFigureOfMeritToCheckSymmetry> lattice_result[NUM_LS]) const
{
ZLOG_INFO( "Determination of the Bravais type is being carried out...\n" );
	putLatticeCandidatesForEachBravaisTypes(S_super, abc_axis, rh_axis, lattice_result);
	for(Int4 i=1; i<NUM_LS; i++)
	{
		for(vector<LatticeFigureOfMeritToCheckSymmetry>::iterator it=lattice_result[i].begin();
				it<lattice_result[i].end(); it++)
		{
			it->setDistance(S_super);
		}
	}
}
