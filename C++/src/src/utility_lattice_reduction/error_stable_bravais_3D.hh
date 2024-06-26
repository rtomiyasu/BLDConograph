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
# ifndef ERROR_STABLE_BRAVAIS_3D_
# define ERROR_STABLE_BRAVAIS_3D_

# include <fstream>
# include <map>
# include "mlist.hh"
# include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
# include "../utility_data_structure/nrutil_nr.hh"
# include "../ControlParam.hh"

using namespace std;

class BravaisLatticeDetermination
{
public:
    class ConstantType
    {
    private:
      static const NRMat<Int4> h_F;
      static const NRMat<Int4> h_I;
      static const NRMat<Int4> h_R;
      static const NRMat<Int4> h_B;

    public:
      static const vector<NRMat<Int4>>& putCmP();
      static const vector<NRMat<Int4>>& putCmC (bool doesPrudentSearch);
      static const vector<NRMat<Int4>>& putCoF();
      static const vector<NRMat<Int4>>& putCoI();
      static const vector<NRMat<Int4>>& putChR (bool doesPrudentSearch);
      
      static const NRMat<Int4>& put_hI () { return h_I; };
      static const NRMat<Int4>& put_hF () { return h_F; };
      static const NRMat<Int4>& put_hR () { return h_R; };
      static const NRMat<Int4>& put_hB () { return h_B; };
  };

  class InputType
  {
  private:
    NRMat<Double> Sobs;
    Double eps;
    bool doesPrudentSearch;
    string axisForBaseCenteredSymmetry;
    string axisForRhombohedralSymmetry;

  public:
    void set(const ControlParam& cData);
    void display();
    inline const NRMat<Double>& put_Sobs() const { return Sobs; };
    inline const Double& put_eps () const { return eps; };
    inline bool put_doesPrudentSearch () const { return doesPrudentSearch; };
    inline const string& put_axisForBaseCenteredSymmetry() const { return axisForBaseCenteredSymmetry; };
    inline const string& put_axisForRhombohedralSymmetry() const { return axisForRhombohedralSymmetry; };
  };

  class ResultType
  {
  public:
	  NRMat<Int4> g;
	  NRMat<Double> S;
	  Double distance;

	  inline void setDistance(const NRMat <Double>& Sobs){ distance = dist(S, mprod (mprod (g, Sobs), transpose (g))); };
	  inline bool operator<(const ResultType& arg) const{ return this->distance < arg.distance; };
  };

private:
    static const vector <string> listNames;
    inline string key2str(const string& key, const string& axisForBaseCenteredSymmetry) const
    {
      static const map <string, string> name2str  {
          {"triclinic"                , "Triclinic"},
          {"primitive monoclinic"     , "Monoclinic(P)"},
          {"base-centered monoclinic" , "Monoclinic"},
          {"primitive orthorhombic"   , "Orthorhombic(P)"},
          {"base-centered orthorhombic" , "Orthorhombic(C)"},
          {"body-centered orthorhombic" , "Orthorhombic(I)"},
          {"face-centered orthorhombic" , "Orthorhombic(F)"},
          {"primitive tetragonal"     , "Tetragonal(P)"},
          {"body-centered tetragonal" , "Tetragonal(I)"},
          {"rhombohedral"             , "Rhombohedral"},
          {"hexagonal"                , "Hexagonal"},
          {"primitive cubic"          , "Cubic(P)"},
          {"body-centered cubic"      , "Cubic(I)"},
          {"face-centered cubic"      , "Cubic(F)"},
      };
      const map <string, string>::const_iterator it=name2str.find (key);
      if ( it == name2str.end() ) return "";
      if ( key == "base-centered monoclinic" )
      {
          string suffix;
          if (axisForBaseCenteredSymmetry == "A") suffix = "(B)";
          else if (axisForBaseCenteredSymmetry == "B") suffix = "(C)";
          else if (axisForBaseCenteredSymmetry == "C") suffix = "(A)";
          return it->second + suffix;
      }
      return it->second;
    };
    InputType m_objInput;
    map<string, vector<ResultType> > m_objResult;

    static vector<ResultType> primitive_monoclinic (const NRMat<Double>& S_obs, const NRMat<Int4>& gb, const Double& eps);
    static vector<ResultType> base_centered_monoclinic (const NRMat<Double>& S_obs, const NRMat<Int4>& gd, const Double& eps,
                    const bool& doesPrudentSearch, const string& axisForBaseCenteredSymmetry);
    static vector<ResultType> face_centered_orthorhombic (const NRMat<Double>& S_obs, const NRMat<Int4>& gd, const Double& eps);
    static vector<ResultType> body_centered_orthorhombic (const NRMat<Double>& S_obs, const NRMat<Int4>& gd2, const Double& eps);
    static vector<ResultType> rhombohedral_centering (const NRMat<Double>& S_obs, const NRMat<Int4>& gd, const Double& eps,
    		        const bool& doesPrudentSearch, const string& axisForRhombohedralSymmetry);
    static vector<ResultType> primitive_monoclinic_to_orthorhombic (const vector<ResultType>& mlist, const Double& eps);
    static vector<ResultType> base_monoclinic_to_orthorhombic (const vector<ResultType>& mlist, const Double& eps,
    				const string& axisForRhombohedralSymmetry);
    static vector<ResultType> orthorhombic_to_tetragonal (const vector<ResultType>& mlist, const Double& eps);
    static vector<ResultType> primitive_monoclinic_to_hexagonal (const vector<ResultType>& ans_mP, const Double& eps);
    static vector<ResultType> orthorhombic_to_cubic (const vector<ResultType>& mlist, const Double& eps);

    void to_output_format(vector<string>& namesOrder, map<string,vector<vector<Double>>>& arg) const;

public:
    static const vector <string>& put_listNames() { return listNames; };

    void set_bravais_class (const InputType& input);
    void reset();
    inline const InputType& putInput() const { return m_objInput; };
    inline const map<string, vector<ResultType> >& putResult() const { return m_objResult; };

    string toText0() const;
    inline void toFile0(const string& savePath) const
    {
    	ofstream outputfile (savePath);
    	outputfile << this->toText0();
    	outputfile.close();
    }

    string toText() const;
    inline void toFile(const string& savePath) const
    {
    	ofstream outputfile (savePath);
    	outputfile << this->toText();
    	outputfile.close();
    }
};

#endif
