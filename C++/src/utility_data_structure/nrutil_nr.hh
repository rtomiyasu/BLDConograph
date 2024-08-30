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
#ifndef _NR_UTIL_H_
#define _NR_UTIL_H_

#include <assert.h>
#include <algorithm>
#include "../utility_func/zstring.hh"

using namespace std;

template <class T>
class NRVec {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRVec();
	explicit NRVec(int n);		// Zero-based array
	NRVec(int n, const T &a);	//initialize to constant value
	NRVec(const NRVec &rhs);	// Copy constructor
	NRVec & operator=(const NRVec &rhs);	//assignment
	NRVec & operator=(const T &a);	//assign a to every element
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	~NRVec();

	NRVec& operator+=(const NRVec& rhs);
	NRVec& operator-=(const NRVec& rhs);
	NRVec& operator*=(const T&);

	NRVec delete_entry(const int& idx) const;
	NRVec insert_entry(const int idx, const T& vec) const;
	string toString() const;
};

template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(0)
{
	assert( n >= 0 );
	if( n > 0 ) v = new T[n];
}

template <class T>
NRVec<T>::NRVec(int n, const T& a) : nn(n), v(0)
{
	assert( n >= 0 );
	if( n > 0 ) v = new T[n];
	for(int i=0; i<n; i++) v[i] = a;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(0)
{
	if( rhs.nn > 0 ) v = new T[nn];
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			delete [] (v);
			v = 0;
			nn=rhs.nn;
			if( nn > 0 ) v= new T[nn];
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}


template <class T>
NRVec<T> & NRVec<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i<nn; i++)
		v[i]=a;
	return *this;
}


template <class T>
inline T & NRVec<T>::operator[](const int i)	//subscripting
{
	assert( 0 <= i && i < nn );
	return v[i];
}

template <class T>
inline const T & NRVec<T>::operator[](const int i) const	//subscripting
{
	assert( 0 <= i && i < nn );
	return v[i];
}

template <class T>
inline int NRVec<T>::size() const
{
	return nn;
}

template <class T>
NRVec<T>::~NRVec()
{
	if (v != 0)
		delete[] (v);
}

template <class T>
NRVec<T>& NRVec<T>::operator+=(const NRVec<T>& rhs)
{
	assert(nn==(int)rhs.size());
	
	for(int k=0; k<nn; k++) v[k]+=rhs[k];
	return *this;
}

template <class T>
NRVec<T> operator+(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	NRVec<T> ans(lhs);
	ans+=rhs;
	return ans;
}

template <class T>
NRVec<T>& NRVec<T>::operator-=(const NRVec<T>& rhs)
{
	assert(nn==(int)rhs.size());
	
	for(int k=0; k<nn; k++) v[k]-=rhs[k];
	return *this;
}

template <class T>
NRVec<T> operator-(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	NRVec<T> ans(lhs);
	ans-=rhs;
	return ans;
}

template <class T>
NRVec<T>& NRVec<T>::operator*=(const T& rhs)
{
	for(int k=0; k<nn; k++) v[k]*=rhs;
	return *this;
}

template <class T>
NRVec<T> operator*(const NRVec<T>& lhs, const T& rhs)
{
	NRVec<T> ans(lhs);
	ans*=rhs;
	return ans;
}

template <class T>
NRVec<T> NRVec<T>::delete_entry(const int& idx) const
{
	const int nsize = this->nn-1;
	NRVec<T> ans (nsize);
	for (int i = 0; i < nsize; i++)
    {
		const int i2=(i>=idx?i+1:i);
		ans[i] = (*this)[i2];
    }
	return ans;
}

template <class T>
NRVec<T> NRVec<T>::insert_entry(const int idx, const T& v) const
{
	assert (0 <= idx && idx <= this->nn);
	const int nsize = this->nn+1;
	NRVec<T> ans (nsize);
	for (int i = 0; i < this->nn; i++)
    {
		const int i2=(i>=idx?i+1:i);
		ans[i2] = (*this)[i];
    }
	ans[idx] = v;
	return ans;
}

template<class T>
string NRVec<T>::toString() const
{
    string ans;
    for (int i=0; i < this->nn; i++)
    {
        ans += num2str((*this)[i], 8) + (i+1>=this->nn?"\n":", ");
    }
    return ans;
}

template <class T>
T inner_product(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	const int isize = lhs.size();
	assert( isize == (int)rhs.size() );
	
	T ans = 0;
	for(int k=0; k<isize; k++) ans += lhs[k]*rhs[k];
	return ans;
}

template <class T>
NRVec<T> product_each_element(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	const int isize = lhs.size();
	assert( isize == (int)rhs.size() );
	
	NRVec<T> ans(isize);
	for(int k=0; k<isize; k++) ans[k] = lhs[k]*rhs[k];
	return ans;
}




template <class T>
class NRMat {
private:
	int nn;
	int mm;
	T **v;
public:
	NRMat();
	NRMat(int n, int m);			// Zero-based array
	NRMat(int n, int m, const T& a);	//Initialize to constant
	NRMat(const NRMat &rhs);		// Copy constructor
	NRMat& operator=(const NRMat& rhs);	//assignment
	NRMat& operator=(const T& a);		//assign a to every element
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	~NRMat();

	NRMat& operator+=(const NRMat&);
	NRMat& operator-=(const NRMat&);
	NRMat& operator*=(const T&);

	NRMat delete_row(const int& idx) const;
	NRMat delete_column(const int& idx) const;
	NRMat insert_row(const int idx, const NRVec<T>& vec) const;
	NRMat insert_column(const int idx, const NRVec<T>& vec) const;
	string toString() const;
};

template <class T>
NRMat<T>::NRMat() : nn(0), mm(0), v(0) {}

template <class T>
NRMat<T>::NRMat(int n, int m) : nn(n), mm(m), v(0)
{
	assert( n >= 0 );
	assert( m >= 0 );
	if(n > 0)
	{
		v = new T*[n];
		if(m > 0)
		{
			v[0] = new T[m*n];
			for (int i=1; i< n; i++)
				v[i] = v[i-1] + m;
		}
		else for (int i=0; i< n; i++) v[i] = 0;
	}
}

template <class T>
NRMat<T>::NRMat(int n, int m, const T &a) : nn(n), mm(m), v(0)
{
	assert( n >= 0 );
	assert( m >= 0 );

	int i,j;
	if(n > 0)
	{
		v = new T*[n];
		if(m > 0)
		{
			v[0] = new T[m*n];
			for (i=1; i< n; i++)
				v[i] = v[i-1] + m;
			for (i=0; i< n; i++)
				for (j=0; j<m; j++)
					v[i][j] = a;
		}
		else for (int i=0; i< n; i++) v[i] = 0;
	}
}

template <class T>
NRMat<T>::NRMat(const NRMat &rhs) : nn(rhs.nn), mm(rhs.mm), v(0)
{
	int i,j;
	if(nn > 0)
	{
		v = new T*[nn];
		if(mm > 0)
		{
			v[0] = new T[mm*nn];
			for (i=1; i< nn; i++)
				v[i] = v[i-1] + mm;
			for (i=0; i< nn; i++)
				for (j=0; j<mm; j++)
					v[i][j] = rhs[i][j];
		}
		else for (int i=0; i<nn; i++) v[i] = 0;
	}
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const NRMat<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j;
		if (nn != rhs.nn || mm != rhs.mm) {
			if(nn > 0) delete[] (v[0]);
			delete[] (v);
			v = 0;
			nn=rhs.nn;
			mm=rhs.mm;
			if( nn > 0 )
			{
				v = new T*[nn];
				if( mm > 0 )
				{
					v[0] = new T[mm*nn];
					for (i=1; i<nn; i++)
						v[i] = v[i-1] + mm;
				}
				else for (i=0; i< nn; i++) v[i] = 0;
			}
		}
		for (i=0; i<nn; i++)
			for (j=0; j<mm; j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}


template <class T>
NRMat<T> & NRMat<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i<nn; i++)
		for (int j=0; j<mm; j++)
			v[i][j] = a;
	return *this;
}

template <class T>
inline T* NRMat<T>::operator[](const int i)	//subscripting: pointer to row i
{
	assert( 0 <= i && i < nn );
	return v[i];
}

template <class T>
inline const T* NRMat<T>::operator[](const int i) const
{
	assert( 0 <= i && i < nn );
	return v[i];
}



template <class T>
inline int NRMat<T>::nrows() const
{
	return nn;
}

template <class T>
inline int NRMat<T>::ncols() const
{
	return mm;
}

template <class T>
NRMat<T>::~NRMat()
{
	if(nn > 0) delete[] (v[0]);
	delete[] (v);
}

template <class T>
NRMat<T>& NRMat<T>::operator+=(const NRMat<T>& rhs)
{
	assert(nn==rhs.nrows());
	assert(mm==rhs.ncols());

	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++) v[k][j] += rhs[k][j];
	return *this;
}

template <class T>
NRMat<T>& NRMat<T>::operator-=(const NRMat<T>& rhs)
{
	assert(nn==rhs.nrows());
	assert(mm==rhs.ncols());

	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++) v[k][j] -= rhs[k][j];
	return *this;
}

/*
template <class T>
NRMat<T> operator-(const NRMat<T>& rhs)
{
	const int nn =rhs.nrows();
	const int mm =rhs.ncols();

	NRMat<T> ans(nn, mm);
	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++) ans[k][j] = -rhs[k][j];
	return ans;
}
*/

template <class T>
NRMat<T>& NRMat<T>::operator*=(const T& rhs)
{
	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++) v[k][j] *= rhs;
	return *this;
}


//template <class T>
//NRMat<T>& NRMat<T>::operator/=(const T& rhs)
//{
//	for(int k=0; k<nn; k++)
//		for(int j=0; j<mm; j++) v[k][j] /= rhs;
//	return *this;
//}


template <class T>
NRMat<T> operator-(const NRMat<T>& lhs, const NRMat<T>& rhs)
{
	NRMat<T> ans = lhs;
	ans -= rhs;
	return ans;
}

template <class T>
NRMat<T> operator*(const NRMat<T>& lhs, const T& rhs)
{
	NRMat<T> ans = lhs;
	ans *= rhs;
	return ans;
}

template<class T>
NRMat<T> operator/(const NRMat<T>& lhs, const int& rhs)
{
	const int nn = lhs.nrows();
	const int mm = lhs.ncols();
	NRMat<T> ans = lhs;
	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++)
		{
			assert( ans[k][j] % rhs == 0 );
			ans[k][j] /= rhs;
		}
	return ans;
}

template <class T>
NRMat<T> NRMat<T>::delete_row(const int& idx) const
{
	const int nrows = this->nn-1;
	const int ncols = this->mm;
	NRMat<T> ans (nrows, ncols);
	for (int i = 0; i < nrows; i++)
    {
		const int i2=(i>=idx?i+1:i);
		for (int j = 0; j < ncols; j++)
	    {
			ans[i][j] = (*this)[i2][j];
	    }
    }
	return ans;
}

template <class T>
NRMat<T> NRMat<T>::delete_column(const int& idx) const
{
	const int nrows = this->nn;
	const int ncols = this->mm-1;
	NRMat<T> ans (nrows, ncols);
	for (int i = 0; i < nrows; i++)
    {
		for (int j = 0; j < idx; j++)
	    {
			ans[i][j] = (*this)[i][j];
	    }
		for (int j = idx, j2 = idx+1; j < ncols; j++, j2++)
	    {
			ans[i][j] = (*this)[i][j2];
	    }
    }
	return ans;
}

template <class T>
NRMat<T> NRMat<T>::insert_row(const int idx, const NRVec<T>& vec) const
{
	assert (0 <= idx && idx <= this->nn);
	assert (this->mm == (int)vec.size());
	const int nrows = this->nn+1;
	NRMat<T> ans (nrows, this->mm);
	for (int i = 0; i < this->nn; i++)
    {
		const int i2=(i>=idx?i+1:i);
		for (int j=0; j < this->mm; j++)
	    {
			ans[i2][j] = (*this)[i][j];
	    }
    }
	for (int j=0; j < this->mm; j++)
    {
		ans[idx][j] = vec[j];
    }
	return ans;
}

template <class T>
NRMat<T> NRMat<T>::insert_column(const int idx, const NRVec<T>& vec) const
{
	assert (0 <= idx && idx <= this->mm);
	assert (this->nn == (int)vec.size());
	const int ncols = this->mm+1;
	NRMat<T> ans (this->nn, ncols);
	for (int i = 0; i < this->nn; i++)
    {
		for (int j=0; j < idx; j++)
	    {
			ans[i][j] = (*this)[i][j];
	    }
		for (int j=idx+1, j2=idx; j < ncols; j++, j2++)
	    {
			ans[i][j] = (*this)[i][j2];
	    }
    }
	for (int i = 0; i < this->nn; i++)
    {
		ans[i][idx] = vec[i];
    }
	return ans;
}

template<class T>
string NRMat<T>::toString() const
{
    string ans;
    for (int i = 0; i < this->nrows(); i++)
    {
        for (int j = 0; j < this->ncols(); j++)
        {
            ans += num2str((*this)[i][j], 8) + (j+1>=this->ncols()?"\n":",");
        }
    }
    return ans;
}


inline NRMat<double> int2double (const NRMat<int> mat)
{
	int irow = mat.nrows();
	int icol = mat.ncols();

	NRMat<double> ans (irow, icol, 0.0);
	for (int i=0; i<irow; i++)
	{
		for (int j=0; j<icol; j++)
		{
			ans[i][j] = double (mat[i][j]);
		}
	}
	return ans;
}

template <class T>
NRMat<T> mprod(const NRMat<T>& lhs, const NRMat<T>& rhs)
{
	const int isize = lhs.ncols();
	assert( isize == rhs.nrows() );

	const int irow = lhs.nrows();
	const int icol = rhs.ncols();
	
	NRMat<T> ans(irow, icol, 0);
	for(int k=0; k<irow; k++)
		for(int j=0; j<icol; j++)
			for(int l=0; l<isize; l++) ans[k][j] += lhs[k][l] * rhs[l][j];
	
	return ans;
}

// ----------over load of mprod in nrutil_nr.cc-------------
// lhs     rhs        output
// int,    int    ---> int
// double, double ---> double
// int,    double ---> double
// double, int    ---> double 
// <<int @ double --> double>>

inline NRMat<double> mprod (const NRMat<int>& lhs, const NRMat<double>& rhs)
{
	NRMat<double> lhs_ = int2double (lhs);
	NRMat<double> ans = mprod (lhs_, rhs);
	return ans;
}

// <<double @ int --> double>>
inline NRMat<double> mprod (const NRMat<double>& lhs, const NRMat<int>& rhs)
{
	NRMat<double> rhs_ = int2double (rhs);
	NRMat<double> ans = mprod (lhs, rhs_);
	return ans;
}

template<class T>
NRMat<T> diagonal (const vector <T> dvec)
{
	int ndim = dvec.size();
	NRMat<T> mat (ndim, ndim, 0);
	for (int i = 0; i < ndim; i++)
	{
		mat[i][i] = dvec[i];
	}
	return mat;
} 

template<class T>
inline T trace (const NRMat<T> mat)
{
	T ans = 0;
	int ndim = mat.nrows();
	for (int i = 0; i < ndim; i++)
	{
		ans += mat[i][i];
	}
	return ans;
}


//template <class T>
//NRMat<T> transpose(const NRMat<T>& lhs)
//{
//	const int icol = lhs.nrows();
//	const int irow = lhs.ncols();
//	
//	NRMat<T> ans(irow, icol);
//	for(int k=0; k<irow; k++)
//		for(int j=0; j<icol; j++) ans[k][j] = lhs[j][k];
//	
//	return ans;
//}

template <class T>
inline NRMat<T> transpose(const NRMat<T>& rhs)
{
	NRMat<T> ans(rhs.ncols(), rhs.nrows());
	for(int i=0; i<rhs.nrows(); i++)
	{
		for(int j=0; j<rhs.ncols(); j++)
		{
			ans[j][i] = rhs[i][j];
		}
	}
	return ans;
}


template <class T>
NRVec<T> left_act(const NRMat<T>& lhs, const NRVec<T>& rhs)
{
	const int irow = lhs.nrows();
	const int icol = lhs.ncols();
	assert(icol==(int)rhs.size());

	
	NRVec<T> ans( irow, 0 );
	for(int k=0; k<irow; k++)
		for(int j=0; j<icol; j++) ans[k] += lhs[k][j]*rhs[j];

	return ans;
}

template <class T>
NRVec<T> right_act(const NRVec<T>& lhs, const NRMat<T>& rhs)
{
	const int irow = rhs.nrows();
	const int icol = rhs.ncols();
	assert(irow==(int)lhs.size());

	
	NRVec<T> ans( icol, 0 );
	for(int k=0; k<icol; k++)
		for(int j=0; j<irow; j++) ans[k] += lhs[j]*rhs[j][k];

	return ans;
}

template<class T>
inline void transpose_square_matrix(NRMat<T>& rhs)
{
	const int isize = rhs.nrows();
	assert( isize == rhs.ncols() );
	
	for(int i=0; i<isize; i++)
		for(int j=0; j<i; j++) swap(rhs[i][j], rhs[j][i]);
}


template<class T>
inline T Determinant3(const NRMat<T>& rhs)
{
	assert( rhs.nrows() == 3 && rhs.ncols() == 3 );
	
	const T det12 = rhs[1][1]*rhs[2][2]-rhs[1][2]*rhs[2][1];
	const T det12_02 = rhs[1][0]*rhs[2][2]-rhs[1][2]*rhs[2][0];
	const T det12_01 = rhs[1][0]*rhs[2][1]-rhs[1][1]*rhs[2][0];
	
	return rhs[0][0]*det12 - rhs[0][1]*det12_02 + rhs[0][2]*det12_01;
}



template<class T>
inline NRMat<double> Inverse3(const NRMat<T>& rhs)
{
	assert( rhs.nrows() == 3 && rhs.ncols() == 3 );
	
	const T det01 = rhs[0][0]*rhs[1][1]-rhs[0][1]*rhs[1][0];
	const T det02 = rhs[0][0]*rhs[2][2]-rhs[0][2]*rhs[2][0];
	const T det12 = rhs[1][1]*rhs[2][2]-rhs[1][2]*rhs[2][1];
	const T det01_02 = rhs[0][0]*rhs[1][2]-rhs[0][2]*rhs[1][0];
	const T det02_01 = rhs[0][0]*rhs[2][1]-rhs[0][1]*rhs[2][0];
	const T det01_12 = rhs[0][1]*rhs[1][2]-rhs[0][2]*rhs[1][1];
	const T det02_12 = rhs[0][1]*rhs[2][2]-rhs[0][2]*rhs[2][1];
	const T det12_02 = rhs[1][0]*rhs[2][2]-rhs[1][2]*rhs[2][0];
	const T det12_01 = rhs[1][0]*rhs[2][1]-rhs[1][1]*rhs[2][0];
	
	const T det = rhs[0][0]*det12 - rhs[0][1]*det12_02 + rhs[0][2]*det12_01;
	assert( det != 0. );
	
	NRMat<T> ans(3,3);
	ans[0][0] = det12;
	ans[0][1] = -det02_12;
	ans[0][2] = det01_12;
	ans[1][0] = -det12_02;
	ans[1][1] = det02;
	ans[1][2] = -det01_02;
	ans[2][0] = det12_01;
	ans[2][1] = -det02_01;
	ans[2][2] = det01;
	return ans * (1./ det);
}

inline NRMat<int> IInverse3(const NRMat<int>& rhs)
{
	assert( rhs.nrows() == 3 && rhs.ncols() == 3 );

	const int det01 = rhs[0][0]*rhs[1][1]-rhs[0][1]*rhs[1][0];
	const int det02 = rhs[0][0]*rhs[2][2]-rhs[0][2]*rhs[2][0];
	const int det12 = rhs[1][1]*rhs[2][2]-rhs[1][2]*rhs[2][1];
	const int det01_02 = rhs[0][0]*rhs[1][2]-rhs[0][2]*rhs[1][0];
	const int det02_01 = rhs[0][0]*rhs[2][1]-rhs[0][1]*rhs[2][0];
	const int det01_12 = rhs[0][1]*rhs[1][2]-rhs[0][2]*rhs[1][1];
	const int det02_12 = rhs[0][1]*rhs[2][2]-rhs[0][2]*rhs[2][1];
	const int det12_02 = rhs[1][0]*rhs[2][2]-rhs[1][2]*rhs[2][0];
	const int det12_01 = rhs[1][0]*rhs[2][1]-rhs[1][1]*rhs[2][0];

	const int det = rhs[0][0]*det12 - rhs[0][1]*det12_02 + rhs[0][2]*det12_01;
	assert( abs(det) == 1 );
	
	NRMat<int> ans(3,3);
	ans[0][0] = det12;
	ans[0][1] = -det02_12;
	ans[0][2] = det01_12;
	ans[1][0] = -det12_02;
	ans[1][1] = det02;
	ans[1][2] = -det01_02;
	ans[2][0] = det12_01;
	ans[2][1] = -det02_01;
	ans[2][2] = det01;

	if( det == 1 )
	{
		return ans;
	}
	else
	{
		return ans * (-1);
	}
}

/*
template <class T>
NRMat<T> transpose_product(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	const int irow = lhs.size();
	const int icol = lhs.size();
	
	NRMat<T> ans( irow, icol );
	for(int k=0; k<irow; k++)
		for(int j=0; j<icol; j++) ans[k][j] = lhs[k]*rhs[j];

	return ans;
}
*/

template <class T>
inline NRMat<T> identity_matrix(const int& isize)
{
	NRMat<T> TransMat(isize,isize,0);
	for(int k=0; k<isize; k++) TransMat[k][k] = 1;

	return TransMat;
}

template <class T>
vector <int> argsort (const vector<T>& vec)
{
    int vecSize = vec.size();
	vector <int> index (vec.size());
	for (int i = 0; i < vecSize; i++) index[i] = i;
    sort (index.begin(), index.end(),
            [&](int x, int y){return vec[x] < vec[y];});
    return index;
}

template <class T>
NRMat<T> put_matrix (const vector<vector<T>>& vec)
{
	int nrow = vec.size();
	int ncol = vec[0].size();
	NRMat<T> ans (nrow, ncol, 0);
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < ncol; j++)
		{
			ans[i][j] = vec[i][j];
		}
	}
	return ans;
}

#endif /* _NR_UTIL_H_ */
