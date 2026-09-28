/*=============================================================================
*   
*   Filename : MyVector.H
*   Creator : Han Zhou
*   Date : 06/15/21
*   Description : 
*
=============================================================================*/
 
#ifndef _MYVECTOR_H
#define _MYVECTOR_H

#include <iostream>
#include <cfloat>
#include <cassert>

template<typename T>
class Vector
{
	protected :

		int _size = 0;
		T *_data = 0;

  protected :
    void allocateMemory(int n)
    {
      _data = new T[n];
    }
    void freeMemory(void)
    {
      delete[] _data;
      _data = 0;
    }

	public :

		Vector(void);
		Vector(int n);
		Vector(int n, const T&);
		Vector(double low, double high, int n);
		Vector(const Vector &v);
		virtual ~Vector(void);

	public :

		T& operator[](int);
		const T &operator[](int) const;

		Vector &operator = (const Vector &v);

		int size(void) const;

		void clear(void);
		void fill(const T &a);
		void allocate(int n);
		void reallocate(int n);
    void assign(int n, const T&a);
		void setAxis(double low, double high, int n);

#ifdef DATA_DEBUG
		void print(void) const;
#endif

		T max(void) const;
		T min(void) const;
		T mean(void) const;
    T sum(void) const;
    T sum0(void) const;
		T max_norm(void) const;
		T l1_norm(void) const;
		T l2_norm(void) const;
		T l1_norm_scaled(void) const;
		T l2_norm_scaled(void) const;

		T max_error(const Vector &u) const;
		T l1_error_scaled(const Vector &u) const;
		T l2_error_scaled(const Vector &u) const;

		T inner_product(const Vector &u) const;

		void normalize(void);

    T* data(void) {return _data;};
    T* getDataPointer(void) {return _data;};
    T* const data(void) const {return _data;};
    T* const getDataPointer(void) const {return _data;};

		/////////////////////////////////////////////////////////////////////

		Vector &operator += (const Vector &v)
		{
			assert(this->_size == v._size);
			for(int d = 0; d < _size; d++){
				this->_data[d] += v._data[d];
			}
			return *this;
		}

		Vector &operator -= (const Vector &v)
		{
			assert(this->_size == v._size);
			for(int d = 0; d < _size; d++){
				this->_data[d] -= v._data[d];
			}
			return *this;
		}

		Vector &operator += (double a)
		{
			for(int d = 0; d < _size; d++){
				this->_data[d] += a;
			}
			return *this;
		}
		Vector &operator -= (double a)
		{
			for(int d = 0; d < _size; d++){
				this->_data[d] -= a;
			}
			return *this;
		}

		Vector &operator += (int a)
		{
			for(int d = 0; d < _size; d++){
				this->_data[d] += a;
			}
			return *this;
		}
		Vector &operator -= (int a)
		{
			for(int d = 0; d < _size; d++){
				this->_data[d] -= a;
			}
			return *this;
		}

		Vector &operator *= (double s)
		{
			for(int d = 0; d < _size; d++){
				this->_data[d] *= s;
			}
			return *this;
		}

		Vector &operator /= (double s)
		{
			assert(fabs(s) >= DBL_MIN);
			for(int d = 0; d < _size; d++){
				this->_data[d] /= s;
			}
			return *this;
		}

		Vector cwiseProduct(const Vector &v) const
		{
			assert(_size == v._size);
			Vector v1(_size);
			for(int d = 0; d < _size; d++){
				v1._data[d] = _data[d] * v._data[d];
			}
			return v1;
		}

    friend double fabs(const Vector &a)
    {
      double res = 0.0;
      for(int d = 0; d < a._size; d++){
        res = std::max(res, fabs(a._data[d]));
      }
      return res;
    }

		friend Vector operator - (const Vector &a, const Vector &b)
		{
      assert(a._size == b._size);

			Vector v1(a._size);
			for(int d = 0; d < a._size; d++){
				v1._data[d] = a._data[d] - b._data[d];
			}
			return v1;
		}

		friend Vector operator + (const Vector &a, const Vector &b)
		{
      assert(a._size == b._size);

			Vector v1(a._size);
			for(int d = 0; d < a._size; d++){
				v1._data[d] = a._data[d] + b._data[d];
			}
			return v1;
		}

		friend T operator * (const Vector &a, const Vector &b)
		{
      assert(a._size == b._size);
      T sum(0);
			for(int d = 0; d < a._size; d++){
        sum += a._data[d]*b._data[d];
			}
			return sum;
		}

		friend Vector operator + (double s, const Vector &a)
		{
			Vector v1(a._size);
			for(int d = 0; d < a._size; d++){
				v1._data[d] = s + a._data[d];
			}
			return v1;
		}

		friend Vector operator * (double s, const Vector &a)
		{
			Vector v1(a._size);
			for(int d = 0; d < a._size; d++){
				v1._data[d] = s * a._data[d];
			}
			return v1;
		}

		friend Vector operator * (const Vector &a, double s)
		{
			Vector v1(a._size);
			for(int d = 0; d < a._size; d++){
				v1._data[d] = s*a._data[d];
			}
			return v1;
		}

		friend std::ostream &operator << (std::ostream &os, const Vector &a)
		{
			for(int i = 0; i < a._size; i++){
				os << a._data[i] << std::endl;
			}
			return os;
		}

		friend T inner_product(const Vector &a, const Vector &b)
    {
      assert(a._size == b._size);
      T sum(0);
			for(int d = 0; d < a._size; d++){
        sum += a._data[d]*b._data[d];
			}
			return sum;
    }

		friend T max_error(const Vector &a, const Vector &b)
		{
      assert(a._size == b._size);

      T temp(0);
			for(int d = 0; d < a._size; d++){
        T err = fabs(a._data[d] - b._data[d]);
        temp = err > temp ? err : temp;
			}
			return temp;
		}

		friend T root_mean_square_error(const Vector &a, const Vector &b)
		{
      assert(a._size == b._size);

      T temp(0);
			for(int d = 0; d < a._size; d++){
        T err = a._data[d] - b._data[d];
        temp += err*err;
			}
      temp = sqrt(temp/a._size);
			return temp;
		}
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(void):_size(0), _data(0)
{
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(int n):_size(n)
{
	_data = new T[_size];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(double low, double high, int n) : _size(n+1)
{
	assert(high > low);
	_data = new T[_size];

	double h = (high - low) / n;
	for(int i = 0; i <= n; i++){
		_data[i] = low + i*h;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(int n, const T &a):_size(n)
{
	_data = new T[_size];
	for(int i = 0; i < n; i++){
		_data[i] = a;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::Vector(const Vector<T> &v):_size(v._size)
{
	allocate(v._size);
	for(int i = 0; i < v._size; i++){
		_data[i] = v._data[i];
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::allocate(int n)
{
	reallocate(n);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::reallocate(int n)
{
	assert(n >= 0);

  if (0 != _data) {
    if (n == _size) {
      return;
    } else {
      freeMemory();
    }
  }
  _size = n;
  allocateMemory(n);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::assign(int n, const T&a)
{
  reallocate(n);
  fill(a);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::setAxis(double low, double high, int n)
{
	assert(high > low);
	reallocate(n+1);

	double h = (high - low) / n;
	for(int i = 0; i <= n; i++){
		_data[i] = low + i*h;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T>::~Vector(void)
{
	clear();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline Vector<T> &Vector<T>::operator = (const Vector<T> &v)
{
	if(this != &v) {
		reallocate(v._size);
		for(int i = 0; i < _size; i++){
			_data[i] = v._data[i];
		}
	}
	return *this;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::clear(void)
{
	if(0 != _data) {
    freeMemory();
	}
	_size = 0;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::fill(const T &a)
{
	for(int i = 0; i < _size; i++){
		_data[i] = a;
	}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline int Vector<T>::size(void) const
{
	return _size;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline const T &Vector<T>::operator[](int i) const
{
	assert(i >= 0 && i < _size);
#ifdef DATA_DEBUG
	if(i < 0 || i >= _size) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}
#endif
	return _data[i];
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T &Vector<T>::operator[](int i)
{
	assert(i >= 0 && i < _size);
#ifdef DATA_DEBUG
	if(i < 0 || i >= _size) {
		std::cout << "Index out of range." << std::endl;
		exit(1);
	}
#endif
	return _data[i];
}

#ifdef DATA_DEBUG
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::print(void) const
{
	if (std::is_same<T, double>::value || 
			std::is_same<T, float>::value || 
			std::is_same<T, long double>::value) {

  	for (int i = 0; i < _size; i++) {
  	  if ((i > 0) && (i%8 == 0)) {
  	    printf("\n"); 
  	  }
  	  printf("%9.6f ", _data[i]); 
  	}
  	printf("\n");

	} else if (std::is_same<T, int>::value || 
						 std::is_same<T, unsigned int>::value || 
						 std::is_same<T, long unsigned int>::value) {

  	for (int i = 0; i < _size; i++) {
  	  if ((i > 0) && (i%16 == 0)) {
  	    printf("\n"); 
  	  }
  	  printf("%4d ", _data[i]); 
  	}
  	printf("\n"); 

	}
}
#endif

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::max(void) const
{
	T max_u = _data[0];
	for(int i = 1; i < _size; i++){
		max_u = max_u > _data[i] ? max_u : _data[i];
	}
	return max_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::min(void) const
{
	T min_u = _data[0];
	for(int i = 1; i < _size; i++){
		min_u = min_u < _data[i] ? min_u : _data[i];
	}
	return min_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::mean(void) const
{
	return (1.0/_size) * sum();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::sum0(void) const
{
  T sum(0); 
  for(int i = 0; i < _size; i++){
    sum += _data[i];
  }
  return sum;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::sum(void) const
{
  T sum(0); 
  T comp(0); 

  for(int i = 0; i < _size; i++){
    T x = _data[i];
    T t = sum + x;
    //comp += (sum - t) + x;
    if (fabs(sum) >= fabs(x)) {
      comp += (sum - t) + x;
    } else {
      comp += (x - t) + sum;
    }
    sum = t;
  }
  return sum + comp;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::max_norm(void) const
{
	T max_u = fabs(_data[0]);
	for(int i = 1; i < _size; i++){
		T u = fabs(_data[i]);
		max_u = max_u > u ? max_u : u;
	}
	return max_u;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l1_norm(void) const
{
	T sum = 0.0;
	for(int i = 0; i < _size; i++){
		sum += fabs(_data[i]);
	}
	return sum;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l2_norm(void) const
{
	T sum = 0.0;
	for(int i = 0; i < _size; i++){
		sum += _data[i] * _data[i];
	}
	return sqrt(sum);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l1_norm_scaled(void) const
{
	T sum = 0.0;
	for(int i = 0; i < _size; i++){
		sum += fabs(_data[i]);
	}
	return sum / _size;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l2_norm_scaled(void) const
{
	T sum = 0.0;
	for(int i = 0; i < _size; i++){
		sum += _data[i] * _data[i];
	}
	return sqrt(sum / _size);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline void Vector<T>::normalize(void)
{
	T sum = 0.0;
	for(int i = 0; i < _size; i++){
		sum += _data[i] * _data[i];
	}

	assert(sum >= 1.0e-15);

#ifdef DATA_DEBUG
	if(sum < 1.0E-15) {
		std::cout << "zero vector." << std::endl;
		exit(1);
	}
#endif

	T r_norm = 1.0 / sqrt(sum);

	for(int i = 0; i < _size; i++){
		_data[i] *= r_norm;
	}

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::max_error(const Vector<T> &u) const
{
	assert(_size == u.size());
#ifdef DATA_DEBUG
	if(_size != u.size()) {
		std::cout << "wrong size. " << std::endl;
		exit(1);
	}
#endif

	T max_e = fabs(_data[0] - u[0]);
	for(int i = 1; i < _size; i++){
		T e = fabs(_data[i] - u[i]);
		max_e = max_e > e ? max_e : e;
	}
	return max_e;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l1_error_scaled(const Vector<T> &u) const
{
	assert(_size == u.size());
#ifdef DATA_DEBUG
	if(_size != u.size()) {
		std::cout << "wrong size. " << std::endl;
		exit(1);
	}
#endif

	T sum = 0.0;
	for(int i = 0; i < _size; i++){
		sum += fabs(_data[i] - u[i]);
	}
	return sum / _size;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename T>
inline T Vector<T>::l2_error_scaled(const Vector<T> &u) const
{
	assert(_size == u.size());
#ifdef DATA_DEBUG
	if(_size != u.size()) {
		std::cout << "wrong size. " << std::endl;
		exit(1);
	}
#endif

	T sum = 0.0;
	for(int i = 0; i < _size; i++){
		T err = _data[i] - u[i];
		sum += err * err;
	}
	return sqrt(sum / _size);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template< typename T >
inline T Vector< T >::inner_product(const Vector<T> &u) const
{
	T sum = 0;
	for(int i = 0; i < _size; i++){
		sum += u[i] * _data[i];
	}
	return sum;
}

//=============================================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif
