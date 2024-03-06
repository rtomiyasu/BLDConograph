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
#include <string>
#include "ControlParam.hh"
#include "utility_func/zstring.hh"
#include "zlog/zlog.hh"


const map<eABCaxis, string> ControlParam::ABCaxisLabel = putABCaxisLabel();
const map<eRHaxis, string> ControlParam::RHaxisLabel = putRHaxisLabel();

const pair< RWParamProperty, RWParamData<void*> > ControlParam::LatticeParameter_Data(
RWParamProperty(VOIDDATA, "LatticeParameter"),
RWParamData<void*>(NULL, REPLACE_NONE<void*>, NULL, NULL, NULL, NULL, -1, -1) );

const pair< RWParamProperty, RWParamData<Double> > ControlParam::LatticeParameters_Data[6]
   	= {
   			pair< RWParamProperty, RWParamData<Double> >(
   					RWParamProperty(DVALUE, "a"),
   					RWParamData<Double>(1.0, REPLACE_NONE<Double>, GT<Double>, 0.0, NULL, MAX_DP(), -1, -1) ),
   		 	pair< RWParamProperty, RWParamData<Double> >(
   		 			RWParamProperty(DVALUE, "b"),
   					RWParamData<Double>(1.0, REPLACE_NONE<Double>, GT<Double>, 0.0, NULL, MAX_DP(), -1, -1) ),
  			pair< RWParamProperty, RWParamData<Double> >(
   					RWParamProperty(DVALUE, "c"),
   					RWParamData<Double>(1.0, REPLACE_NONE<Double>, GT<Double>, 0.0, NULL, MAX_DP(), -1, -1) ),

  			pair< RWParamProperty, RWParamData<Double> >(
   					RWParamProperty(DVALUE, "alpha"),
   					RWParamData<Double>(90.0, REPLACE_NONE<Double>, GT<Double>, 0.0, LT<Double>, 180.0, -1, -1) ),
  			pair< RWParamProperty, RWParamData<Double> >(
   					RWParamProperty(DVALUE, "beta"),
   					RWParamData<Double>(90.0, REPLACE_NONE<Double>, GT<Double>, 0.0, LT<Double>, 180.0, -1, -1) ),
  			pair< RWParamProperty, RWParamData<Double> >(
   					RWParamProperty(DVALUE, "gamma"),
   					RWParamData<Double>(90.0, REPLACE_NONE<Double>, GT<Double>, 0.0, LT<Double>, 180.0, -1, -1) )
	};

const pair<RWParamProperty, RWParamData<bool> > ControlParam::DoesPrudentSymSearch_Data(
		RWParamProperty(BOOLFLAG, "DoesPrudentSymmetrySearch"),
		RWParamData<bool>(false, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) );
const pair<RWParamProperty, RWParamData<bool> > ControlParam::OutputSymmetry_Data[NUM_LS]
	= {
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputTriclinic"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputMonoclinicP"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputMonoclinicB"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputOrthorhombicP"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputOrthorhombicB"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputOrthorhombicI"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputOrthorhombicF"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputTetragonalP"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputTetragonalI"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputRhombohedral"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputHexagonal"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputCubicP"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputCubicI"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputCubicF"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) )
		};

const pair<RWParamProperty, RWParamData<Double> > ControlParam::Resol_Data(
		RWParamProperty(DVALUE, "Resolution"), 
		RWParamData<Double>(0.05, REPLACE_NONE<Double>, GE<Double>, 0.0, LE<Double>, 0.25, -1, -1) );	// 0 <= param < INF.

const pair<RWParamProperty, RWParamData<string> > ControlParam::MonoBaseAxis_Data(
		RWParamProperty(STRVALUE, "AxisForBaseCenteredSymmetry"),
		RWParamData<string>("B", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );
const pair<RWParamProperty, RWParamData<string> > ControlParam::RhomAxis_Data(
		RWParamProperty(STRVALUE, "AxisForRhombohedralSymmetry"),
		RWParamData<string>("Hexagonal", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );

const Int4 ControlParam::NumCores = 1;

ControlParam::ControlParam()
	: m_length(LatticeParameters_Data[0].second.initial_value, LatticeParameters_Data[1].second.initial_value, LatticeParameters_Data[2].second.initial_value),
     	m_angle(LatticeParameters_Data[3].second.initial_value, LatticeParameters_Data[4].second.initial_value, LatticeParameters_Data[5].second.initial_value),
		DoesPrudentSymSearch(DoesPrudentSymSearch_Data.second.initial_value),
		Resol(Resol_Data.second.initial_value),
		MonoBaseAxis(MonoBaseAxis_Data.second.initial_value),
		RhomAxis(RhomAxis_Data.second.initial_value)
{
	for(Int4 i=0; i<NUM_LS; i++)
	{
		OutputSymmetry[i] = OutputSymmetry_Data[i].second.initial_value;
	}
}


ControlParam::~ControlParam()
{
}


const string& ControlParam::putTagLabel() const
{
	static const string label = "BLDConographParameters";
	return label;
}




void ControlParam::setData(const RWParamProperty& parent_prop,
		vector<RWParam_void>& tray)
{
	if( parent_prop.putLabel() == this->putTagLabel() )
	{
		tray.push_back( RWParam_void(LatticeParameter_Data) );

		tray.push_back( RWParam_void(DoesPrudentSymSearch_Data, &DoesPrudentSymSearch) );
		for(Int4 i=0; i<NUM_LS; i++)
		{
			tray.push_back( RWParam_void(OutputSymmetry_Data[i], &OutputSymmetry[i]) );
		}
		tray.push_back( RWParam_void(Resol_Data, &Resol) );
		tray.push_back( RWParam_void(MonoBaseAxis_Data, &MonoBaseAxis) );
		tray.push_back( RWParam_void(RhomAxis_Data, &RhomAxis) );
	}
	else if( IsEqualTag(parent_prop, LatticeParameter_Data.first) )
	{
		tray.push_back( RWParam_void( LatticeParameters_Data[0], &m_length[0] ) );
		tray.push_back( RWParam_void( LatticeParameters_Data[1], &m_length[1] ) );
		tray.push_back( RWParam_void( LatticeParameters_Data[2], &m_length[2] ) );
		tray.push_back( RWParam_void( LatticeParameters_Data[3], &m_angle[0] ) );
		tray.push_back( RWParam_void( LatticeParameters_Data[4], &m_angle[1] ) );
		tray.push_back( RWParam_void( LatticeParameters_Data[5], &m_angle[2] ) );
	}
}


ZErrorMessage ControlParam::checkData(const RWParam_void& param) const
{
	const string Label = param.putLabel(); 
	if( IsEqualTag(param.putProperty(), LatticeParameter_Data.first) )
	{
		SymMat<Double> S(3);
		calCoParameter(m_length, m_angle, S);

ZLOG_INFO( "Input metric tensor S (input) := \n"
			+ num2str( S(0,0) ) + "\n"
			+ num2str( S(1,0) ) + " " + num2str( S(1,1) ) + "\n"
			+ num2str( S(2,0) ) + " " + num2str( S(2,1) ) + " " + num2str( S(2,2) ) + "\n\n" );

		if( Determinant3(S) <= 0.0 )
		{
			return nerror_arg("The metric tensor is not positive definite", __FILE__, __LINE__, __FUNCTION__ );
		}
	}
	else if( IsEqualTag(param.putProperty(), MonoBaseAxis_Data.first) )
	{
		if( find_key(ABCaxisLabel, MonoBaseAxis) == eABCaxis(-1) )
		{
			return nerror_out_range(param.putLabel(), MonoBaseAxis, __FILE__, __LINE__, __FUNCTION__ );
		}
		else return ZErrorMessage();
	}
	else if( IsEqualTag(param.putProperty(), RhomAxis_Data.first) )
	{
		if( find_key(RHaxisLabel, RhomAxis) == eRHaxis(-1) )
		{
			return nerror_out_range(param.putLabel(), RhomAxis, __FILE__, __LINE__, __FUNCTION__ );
		}
		else return ZErrorMessage();
	}
	return I_ReadData::checkData(param);
}
