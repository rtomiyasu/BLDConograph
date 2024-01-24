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
#include <fstream>
#include "p_out.hh"
#include "ControlParam.hh"
#include "lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.hh"
#include "utility_lattice_reduction/matrix_NbyN.hh"


void printHKLdata(const vector<LatticeFigureOfMeritToCheckSymmetry> lattice_result[],
		const ControlParam& cdata,
		const Int4& cal_time,
		const string& fname)
{
	static const Int4 NUM_LS = 14;
	static const string CS_LABEL[NUM_LS] =
			{	"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14" };


	ofstream ofs(fname.c_str());
    ostream *os;
    os = &ofs;
    
   	os->setf(ios::scientific);
   	os->precision(4);
   	
	pair< vector<Int4>::const_iterator, vector<Int4>::const_iterator> it_pair;
	VecDat3<Double> length_axis, angle_axis; 

	Int4 label_start=0;
	os->width(label_start);
	*os << "" << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";

	os->width(label_start);
	*os << "" << "<ZCodeParameters>\n";
	label_start++;

	os->width(label_start);
	*os << "" << "<BLDConographOutput>\n";
	label_start++;

	os->width(label_start);
	*os << "" << "<!-- Bravais lattice with minimal distance.\n";
	os->width(label_start);
	*os << "" << "     (The distance between two unit cells is computed by\n";
	os->width(label_start);
	*os << "" << "           Trace((S-T)^2) := \\sum_{1 <= i,j <= 3} (sij-tij)^2 \n";
	os->width(label_start);
	*os << "" << "      using the Delone-reduced metric tensors S := (sij) and T := (tij) of their reciprocal lattices.)\n\n";
	os->width(3);
	*os << "" << "CentringType  NumberOfSolusions  DistanceFromInputUnitCell" << " : unit-cell parameters.\n";
	for(Int4 k=NUM_LS-1; k>=1; k--)
	{
		os->width(18);
		*os << put_bravais_type_name(eBravaisType(k), cdata.putBaseCenteredAxis());
		os->width(5);
		*os << lattice_result[k].size();

		if( lattice_result[k].empty() )
		{
			*os << "\n";
			continue;
		}

		const LatticeFigureOfMerit& lat_fom = lattice_result[k].begin()->putLatticeFigureOfMerit();

	  	os->width(14);
		*os << sqrt( lattice_result[k].begin()->putDistance() );

		*os << " : ";
		lat_fom.putReducedLatticeConstantsDegree(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), length_axis, angle_axis);
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
	  	*os << endl;
	}
	*os << endl;
	os->width(label_start);
	*os << "" << "-->\n\n";

	for(Int4 k=NUM_LS-1; k>=1; k--)
	{
		*os << "\n";
		os->width(label_start);
		*os << "" << "<!-- Candidates for " << put_bravais_type_name(eBravaisType(k), cdata.putBaseCenteredAxis()) << " -->\n\n";
		
		const Int4 num_topo = lattice_result[k].size();

		for(Int4 n=0; n<num_topo; n++)
	   	{
			const LatticeFigureOfMeritToCheckSymmetry& ans = lattice_result[k][n];

			os->width(label_start);
	      	*os << "" << "<UnitCellCandidate number=\"" << CS_LABEL[(size_t)ans.putLatticeFigureOfMerit().enumBravaisType()] + "0" + num2str<Int4>(num_topo+1) << "\">\n";
	      	label_start++;
	      	ans.printLatticeInformation(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), label_start, os);
	      	
	      	label_start--;
			os->width(label_start);
	      	*os << "" << "</UnitCellCandidate>\n\n";
	   	}
   	   	*os << endl;
	}
	
	label_start--;
	os->width(label_start);
	*os << "" << "</BLDConographOutput>\n";

	label_start--;
	os->width(label_start);
	*os << "" << "</ZCodeParameters>\n";
}
