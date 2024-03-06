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
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <time.h>
#include "ControlFile.hh"
#include "ControlParam.hh"
#include "bravais_type/BravaisLattice.hh"
#include "lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.hh"
#include "zerror_type/error_out.hh"
#include "zerror_type/error_mes.hh"
#include "p_out.hh"
#include "utility_func/zstring.hh"
#include "zlog/zlog.hh"


using namespace std;


int main(int argc, char* argv[])
{
	clock_t start;
	start = clock();    /* Record the starting time. */

	static const string controlFile = "cntl.inp.xml";
	static const string InputFileLabel = "ZCodeParameters";
    
	try{
	    CRLog::append(new CCoutListner());
	    CRLog::append(new FileoutListner("LOG_BLDCONOGRAPH.txt", zListnerID(1)));

ZLOG_INFO( "Reading " + controlFile + "...\n\n" );

		ControlFile cf;
    	ZErrorMessage zerr;
    	zerr = cf.readFile(controlFile, InputFileLabel);
    	if( zerr.putErrorType() != ZErrorNoError) throw zerr;


    	ControlParam cData;
    	zerr = cData.readFile(cf.putControlParamFileName(), InputFileLabel);
    	if( zerr.putErrorType() != ZErrorNoError) throw zerr;

		string fname0;
    	removeFileExtension(cf.putOutputFileName(), fname0);

    	SymMat<Double> S_super(4);
    	if( !cData.putDeloneReducedLatticeParameter(S_super) )
    	{
    		throw nerror_arg("Delone reduction was failed for input unit-cell parameters", __FILE__, __LINE__, __FUNCTION__ );
    	}
ZLOG_INFO( "Delone reduced metric tensor Sdel (input) := \n"
    	+ num2str( S_super(0,0) ) + "\n"
		+ num2str( S_super(1,0) ) + " " + num2str( S_super(1,1) ) + "\n"
		+ num2str( S_super(2,0) ) + " " + num2str( S_super(2,1) ) + " " + num2str( S_super(2,2) ) + "\n\n" );

#ifdef _OPENMP
		omp_set_num_threads(min(omp_get_max_threads(), cData.putNumberOfThreadsToUse()));
ZLOG_INFO( "The number of threads is set to " + num2str( min(omp_get_max_threads(), cData.putNumberOfThreadsToUse()) ) + "\n" );
#endif

		// Indexing.
		static const Int4 NUM_LS = put_number_of_bravais_types();
		vector<LatticeFigureOfMeritToCheckSymmetry> lattice_result[NUM_LS];
		
		BravaisLattice srl;
		srl.setParam(cData);

		start = clock();    /* Record the starting time. */
		srl.execute(S_super, cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), lattice_result);

		// Sort solutions by their distances from the input lattice.
		for(Int4 i=0; i<NUM_LS; i++)
		{
			stable_sort( lattice_result[i].begin(), lattice_result[i].end() );
		}

		// Sets output_flag.
ZLOG_INFO( "Outputting results...\n" );
		const Int4 calculation_time = (clock() - start) / CLOCKS_PER_SEC;

		// Solution having the top M is output as the best solution.
		printHKLdata(lattice_result, cData,
  						calculation_time,
						cf.putOutputFileName());

ZLOG_INFO( "The program has finished Bravais lattice determination in CPU time : " + num2str( Double(clock() - start) / CLOCKS_PER_SEC )
				+ " [sec.]\n" );
	}
	catch(bad_alloc& ball){
		ZErrorMessage zerr = nerror(ball, __FILE__, __LINE__, __FUNCTION__);
		ZLOG_ERROR( zerr.printErrorLog() );
		return 0;
	}
	catch(out_of_range&)
	{
		ZErrorMessage zerr("out_of_range exception has occurred", __FILE__, __LINE__, __FUNCTION__);
		ZLOG_ERROR( zerr.printErrorLog() );
		return 0;
	}
	catch(const ZErrorMessage& etype)
	{
		ZLOG_ERROR( etype.printErrorLog() );
		return 0;
	}

    return 1;
}
