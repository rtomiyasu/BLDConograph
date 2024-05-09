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

#include "ControlFile.hh"
#include "ControlParam.hh"
#include "../src/utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "../src/utility_lattice_reduction/error_stable_bravais_3D.hh"
#include "../src/utility_data_structure/nrutil_nr.hh"

using namespace std;


int main(int argc, char* argv[])
{
	static const string controlFile = "cntl.inp.xml";
	static const string InputFileLabel = "ZCodeParameters";
    
	ControlFile cf;
    ZErrorMessage zerr;
    zerr = cf.readFile(controlFile, InputFileLabel);
    if( zerr.putErrorType() != ZErrorNoError) throw zerr;

    ControlParam cData;
    zerr = cData.readFile(cf.putControlParamFileName(), InputFileLabel);
    if( zerr.putErrorType() != ZErrorNoError) throw zerr;

	cout << "Start processing " + cf.putControlParamFileName() + "..." << endl;
		
	// Set input parameters.
	BravaisLatticeDetermination::InputType bld_input;
	bld_input.set(cData);
	bld_input.display();

	// Run.
	BravaisLatticeDetermination bld;
	bld.set_bravais_class(bld_input);

	cout << "Outputting " + cf.putOutputFileName() + "\n" << endl;
	bld.toFile(cf.putOutputFileName());
	cout << "Completed.\n" << endl;
    return 1;
}
