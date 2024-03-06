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
#include <set>
#include "check_equiv.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "../bravais_type/BravaisType.hh"
#include "lattice_symmetry.hh"

inline bool checkIfFirstEntryIsPositive(const NRMat<Int4>& rhs)
{
	assert( rhs.nrows() == 4 && rhs.ncols() == 3 );
	for(Int4 i=0; i<3; i++)
		for(Int4 j=0; j<3; j++)
		{
			if( rhs[i][j] > 0 ) return true;
			else if( rhs[i][j] < 0 ) return false;
		}
	return false;
}

inline bool operator<(const NRMat<Int4>& lhs, const NRMat<Int4>& rhs)
{
	assert( lhs.nrows() >= 3 && lhs.ncols() == 3 );
	assert( rhs.nrows() == lhs.nrows() && rhs.ncols() == lhs.ncols() );
	for(Int4 i=0; i<3; i++)
		for(Int4 j=0; j<3; j++)
		{
			if( lhs[i][j] < rhs[i][j] ) return true;
			if( rhs[i][j] < lhs[i][j] ) return false;
		}
	return false;
}


static void insert_super_obtuse_equiv(const pair< NRMat<Int4>, SymMat<Double> >& S_super,
		const Double& resol, const Double& Max_totalQ,
		map<NRMat<Int4>, SymMat<Double> >& equivalent_tray)
{
	const Double totalQ = S_super.second(0,0) + S_super.second(1,1) + S_super.second(2,2) + S_super.second(3,3);

	map<NRMat<Int4>, SymMat<Double> >::iterator it;
	Int4 k,l;
	for(Int4 i=0; i<3; i++)
	{
		for(Int4 j=i+1; j<4; j++)
		{
			const Double& Sij = S_super.second(i,j);
			
			if( 0.0 < Sij ) continue; 
			if( !equiv_zero(S_super.second, i,j, resol) ) continue;
			if( Max_totalQ < totalQ - Sij * 2.0 ) continue;
			
			const NRMat<Int4>& mat = put_reduction_matrix_dim_4(i,j);
			pair< NRMat<Int4>, SymMat<Double> > S_super_equiv( mprod(mat, S_super.first),
														transform_sym_matrix(mat, S_super.second) );

			// Ki, Kj, Kk, Kl -> -Ki, Kj, Ki+Kk, Ki+Kl => S_super_equiv.second(k,l) = S_super.second(k,l) - S_super.second(i,j)
			// sum_{i=0}^4 S_super_equiv.second(i,i) = -2*S_super(i,j) + sum_{i=0}^4 S_super.second(i,i)
			put_complement_set4(i,j,k,l);
			if( 0.0 < S_super_equiv.second(k,l) && !equiv_zero(S_super_equiv.second, k,l, resol) ) continue;
			moveSmallerDiagonalLeftUpper(S_super_equiv.second, S_super_equiv.first);

			if( !checkIfFirstEntryIsPositive( S_super_equiv.first ) ) S_super_equiv.first *= -1;

			it = equivalent_tray.find(S_super_equiv.first);
			
			if( it == equivalent_tray.end() )
			{
				equivalent_tray.insert( S_super_equiv );
				insert_super_obtuse_equiv( S_super_equiv, resol, Max_totalQ, equivalent_tray );
			}
		}
	}
};



static void set_super_obtuse_equiv(const SymMat<Double>& S_super_obtuse, const Double& resol,
		map< NRMat<Int4>, SymMat<Double> >& equivalent_tray)
{
	assert( S_super_obtuse.size() == 4 );
	equivalent_tray.clear();
	
	const pair< NRMat<Int4>, SymMat<Double> > S_super_obtuse_pair(put_transform_matrix_row3to4(), S_super_obtuse);
	equivalent_tray.insert( S_super_obtuse_pair );

	const Double Max_totalQ = ( S_super_obtuse(0,0) + S_super_obtuse(1,1)
									+ S_super_obtuse(2,2) + S_super_obtuse(3,3) ) * ( 1.0 + resol );

	insert_super_obtuse_equiv( S_super_obtuse_pair, resol, Max_totalQ, equivalent_tray );
};



inline void permute_row_column(const Int4& i, const Int4& j, pair< NRMat<Int4>, SymMat<Double> >& rhs)
{
	const NRMat<Int4> tmat = put_permutation_matrix_dim_4(i, j);
	rhs.first = mprod( tmat, rhs.first );
	if( !checkIfFirstEntryIsPositive( rhs.first) ) rhs.first *= -1;
	rhs.second = transform_sym_matrix( tmat, rhs.second );
};


// On input, for any i < j such that stage < i, j < 4, S_super(i, i) <= S_super(j, j) 
static void insert_super_permuted(const Int4& stage, 
		const pair< NRMat<Int4>, SymMat<Double> >& S_super,
		const Double& resol,
		map< NRMat<Int4>,  SymMat<Double> >& equivalent_tray)
{
	if( stage >= 3 )
	{
		equivalent_tray.insert( S_super );
		return;
	}
	
	insert_super_permuted(stage+1, S_super, resol, equivalent_tray);

	pair< NRMat<Int4>, SymMat<Double> > S_super_permuted = S_super;
	
	for(Int4 i=stage+1; i<4; i++)
	{
		if( !equiv_resol(S_super.second(stage, stage), S_super.second(i,i), resol ) ) break;
		permute_row_column( stage, i, S_super_permuted );
		insert_super_permuted(stage+1, S_super_permuted, resol, equivalent_tray);
	}
};


static void set_super_permuted(const pair< NRMat<Int4>, SymMat<Double> >& S_super,
		const Double& resol,
		map< NRMat<Int4>,  SymMat<Double> >& equivalent_tray)
{
	equivalent_tray.clear();
	insert_super_permuted( 0, S_super, resol, equivalent_tray );
};



void put_S_super_obtuse_equiv(const SymMat<Double>& S_super_obtuse, const Double& resol,
		vector< SymMat<Double> >& S_super_obtuse_equiv)
{
	S_super_obtuse_equiv.clear();
	
	map< NRMat<Int4>,  SymMat<Double> > super_equiv_tray;
	set_super_obtuse_equiv( S_super_obtuse, resol, super_equiv_tray );

	map< NRMat<Int4>,  SymMat<Double> > permute_equiv_tray;
	map< NRMat<Int4>,  SymMat<Double> >::const_iterator it2;
	map< NRMat<Int4>,  SymMat<Double> > S_super_obtuse_equiv_tray;
	for(map< NRMat<Int4>,  SymMat<Double> >::iterator it=super_equiv_tray.begin(); it!=super_equiv_tray.end(); it++)
	{
		it2 = S_super_obtuse_equiv_tray.find(it->first);
		if( it2 != S_super_obtuse_equiv_tray.end() ) continue;

		set_super_permuted( *it, resol, permute_equiv_tray );
		for(map< NRMat<Int4>,  SymMat<Double> >::const_iterator it3=permute_equiv_tray.begin(); it3!=permute_equiv_tray.end(); it3++)
		{
			S_super_obtuse_equiv_tray.insert(*it3);
			S_super_obtuse_equiv.push_back(it3->second);
		}
	}
}
