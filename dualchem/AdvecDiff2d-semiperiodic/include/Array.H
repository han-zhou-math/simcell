/*=============================================================================
*   
*   Filename : Array.H
*   Creator : Han Zhou
*   Date : 06/15/21
*   Description : 
*
=============================================================================*/
 
#ifndef _ARRAY_H
#define _ARRAY_H

#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

template<typename T, int N>
class Array
{
	protected :

		T _data[N];

	public :

		virtual ~Array(void){};

		Array(void)
		{
			for(int d = 0; d < N; d++){
				_data[d] = 0;
			}
		}

		Array(T a)
		{
			for(int d = 0; d < N; d++){
				_data[d] = a;
			}
		}

		//Array(const T a[N])
		//{
		//	for(int d = 0; d < N; d++){
		//		_data[d] = a[d];
		//	}
		//}

		Array(T a, T b, T c=0, T d=0, T e=0, T f=0)
		{
			T crd[6] = {a, b, c, d, e, f};
			int len = N < 6 ? N : 6;
			for(int m = 0; m < len; m++){
				_data[m] = crd[m];
			}
		}

		Array(const Array& a)
		{
			for(int d = 0; d < N; d++){
				_data[d] = a._data[d];
			}
		}

	public :

		T& operator[](int);
		const T& operator[](int) const;

		Array &operator = (const Array &v)
		{
			if(this != &v) {
				for(int i = 0; i < N; i++){
					_data[i] = v._data[i];
				}
			}
			return *this;
		}

		//Array &operator = (T a)
		//{
		//	for(int i = 0; i < N; i++){
		//		_data[i] = a;
		//	}
		//	return *this;
		//}

		Array operator + (const Array &v)
		{
			Array v1;
			for(int d = 0; d < N; d++){
				v1._data[d] = this->_data[d] + v._data[d];
			}
			return v1;
		}

		Array operator - (const Array &v)
		{
			Array v1;
			for(int d = 0; d < N; d++){
				v1._data[d] = this->_data[d] - v._data[d];
			}
			return v1;
		}

		Array &operator += (const Array &v)
		{
			for(int d = 0; d < N; d++){
				this->_data[d] += v._data[d];
			}
			return *this;
		}

		Array &operator -= (const Array &v)
		{
			for(int d = 0; d < N; d++){
				this->_data[d] -= v._data[d];
			}
			return *this;
		}

		Array &operator += (double a)
		{
			for(int d = 0; d < N; d++){
				this->_data[d] += a;
			}
			return *this;
		}

		Array &operator -= (double a)
		{
			for(int d = 0; d < N; d++){
				this->_data[d] -= a;
			}
			return *this;
		}

		Array &operator *= (double s)
		{
			for(int d = 0; d < N; d++){
				this->_data[d] *= s;
			}
			return *this;
		}

		Array &operator /= (double s)
		{
			assert(fabs(s) >= DBL_MIN);
			for(int d = 0; d < N; d++){
				this->_data[d] /= s;
			}
			return *this;
		}

		Array cwiseProduct(const Array &v) const
		{
			Array v1;
			for(int d = 0; d < N; d++){
				v1._data[d] = _data[d] * v._data[d];
			}
			return v1;
		}

    friend double fabs(const Array &a)
    {
      double res = 0.0;
      for(int d = 0; d < N; d++){
        res = std::max(res, fabs(a._data[d]));
      }
      return res;
    }

		friend Array operator - (const Array &a, const Array &b)
		{
			Array v1;
			for(int d = 0; d < N; d++){
				v1._data[d] = a._data[d] - b._data[d];
			}
			return v1;
		}

		friend Array operator + (const Array &a, const Array &b)
		{
			Array v1;
			for(int d = 0; d < N; d++){
				v1._data[d] = a._data[d] + b._data[d];
			}
			return v1;
		}

		friend T operator * (const Array &a, const Array &b)
		{
      T sum(0);
			for(int d = 0; d < N; d++){
        sum += a._data[d]*b._data[d];
			}
      return sum;
		}

		friend Array operator + (double s, const Array &a)
		{
			Array v1;
			for(int d = 0; d < N; d++){
				v1._data[d] = s + a._data[d];
			}
			return v1;
		}


		friend Array operator * (double s, const Array &a)
		{
			Array v1;
			for(int d = 0; d < N; d++){
				v1._data[d] = s*a._data[d];
			}
			return v1;
		}

		friend Array operator * (const Array &a, double s)
		{
			Array v1;
			for(int d = 0; d < N; d++){
				v1._data[d] = s*a._data[d];
			}
			return v1;
		}

		friend std::ostream &operator << (std::ostream &os, const Array &a)
		{
			for(int d = 0; d < N-1; d++){
				os << a._data[d] << " ";
			}
			os << a._data[N-1];
			return os;
		}

    friend Array normalize(const Array &a)
    {
      auto b = a;
      b.normalize();
      return b;
    }

		void fill(T a)
		{
			for(int d = 0; d < N; d++){
				_data[d] = a;
			}
		}

		void print(void) const;
		void normalize(void);
		T inner_product(const Array &u) const;

		double l2_norm(void) const
		{
			double sum = 0.0;
			for(int d = 0; d < N; d++){
				sum += _data[d]*_data[d];
			}
			return sqrt(sum);
		}

};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, int N>
inline const T& Array<T, N>::operator[](int i) const
{
	assert(i >= 0 && i < N);
	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, int N>
inline T& Array<T, N>::operator[](int i)
{
	assert(i >= 0 && i < N);
	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, int N>
inline void Array<T, N>::normalize(void)
{
	T sum = 0;
	for(int i = 0; i < N; i++){
		sum += _data[i] * _data[i];
	}

	if(sum < 1.0E-15) {
		return;
		//std::cout << "zero array." << std::endl;
		//exit(1);
	}

	T r_n = 1.0 / sqrt(sum);
	for(int i = 0; i < N; i++){
		_data[i] *= r_n;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, int N>
inline T Array<T, N>::inner_product(const Array<T, N> &u) const
{
	T sum = 0;
	for(int i = 0; i < N; i++){
		sum += _data[i] * u[i];
	}
	return sum;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T, int N>
inline void Array<T, N>::print(void) const
{
	std::cout << "(";
	for(int i = 0; i < N - 1; i++){
		std::cout << _data[i] << ", ";
	}
	std::cout << _data[N - 1] << ")" << std::endl;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

inline void CrossProduct(Array<double, 3> const &u, 
                         Array<double, 3> const &v, 
                         Array<double, 3> &w)
{
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - u[0] * v[2];
  w[2] = u[0] * v[1] - u[1] * v[0];
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

