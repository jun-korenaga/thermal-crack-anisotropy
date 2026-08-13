/*
 * array.h
 */

#ifndef _JK_NEW_ARRAY_H_
#define _JK_NEW_ARRAY_H_

#include <vector>
#include <cstddef>
#include <iostream>
#include "emt/error.hpp"

using std::ostream;
using std::size_t;
using std::vector;

template<class T>
class Array1d : public vector<T> {
public:
    Array1d(size_t n=0) : vector<T>(n) {}
    Array1d(const T*, size_t);	// array initialization
    
    T& operator()(size_t i) {
//	if (i>size()){cerr << "Array1d OoR " << i << " " << size() << '\n';}
	return operator[](i); }
    const T& operator()(size_t i) const {
//	if (i>size()){cerr << "Array1d OoR " << i << " " << size() << '\n';}
	return operator[](i); }

    size_t size() const { return vector<T>::size(); }
    T* begin() { return &(vector<T>::front()); }
    T* end() { return (&(vector<T>::back()))+1; }
    const T* begin() const { return &(vector<T>::front()); }
    const T* end() const { return (&(vector<T>::back()))+1; }


    void resize(size_t n, const T& val = T()) { vector<T>::resize(n, val); }
    void push_back (const T& x) { vector<T>::push_back(x); }

    // unary operators
    Array1d<T> operator-();
    
    // binary operators
    Array1d<T>& operator=(const T& val);
    Array1d<T>& operator+=(const Array1d<T>&);
    Array1d<T>& operator+=(const T& val);
    Array1d<T>& operator-=(const Array1d<T>&);
    Array1d<T>& operator-=(const T& val);
    Array1d<T>& operator*=(const T&);
    Array1d<T>& operator/=(const T&);
    
protected:
          T& operator[](size_t i) { return vector<T>::operator[](i-1); }
    const T& operator[](size_t i) const { return vector<T>::operator[](i-1); }
};

template<class T>
inline
Array1d<T>::Array1d(const T* init, size_t n)
    : vector<T>(n)
{
    T* dest = begin()+n;
    const T* src = init+n;
    while (dest > begin()) *--dest = *--src;
}
    
template<class T>
inline
Array1d<T> Array1d<T>::operator-()
{
    Array1d<T> negative(size());

    T* dest = negative.end();
    const T* src = end();
    while (dest > negative.begin()) *--dest = (-1)*(*--src);
    return negative;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest = val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator+=(const Array1d<T>& a)
{
    if (size() != a.size()) error("Array1d::operator+= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest += *--src;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator+=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest += val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator-=(const Array1d<T>& a)
{
    if (size() != a.size()) error("Array1d::operator-= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest -= *--src;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator-=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest -= val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator*=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest *= a;
    return *this;
}
    
template<class T>
inline
Array1d<T>& Array1d<T>::operator/=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest /= a;
    return *this;
}


// nonmember functions
template<class T>
inline
Array1d<T> operator+(const Array1d<T>& a, const Array1d<T>& b)
{
    Array1d<T> r=a;
    return r+=b;
}

template<class T>
inline
Array1d<T> operator-(const Array1d<T>& a, const Array1d<T>& b)
{
    Array1d<T> r=a;
    return r-=b;
}

template<class T>
inline
Array1d<T> operator*(const Array1d<T>& a, const T& val)
{
    Array1d<T> r=a;
    return r*=val;
}

template<class T>
inline
Array1d<T> operator/(const Array1d<T>& a, const T& val)
{
    Array1d<T> r=a;
    return r/=val;
}

template<class T>
ostream& operator<<(ostream& s, const Array1d<T>& a)
{
    for (int i=1; i<a.size(); i++){
	s << a(i) << ", ";
    }
    s << a(a.size()) << '\n';
    return s;
}

// Array2d
template<class T>
class Array2d : public vector<T> {
public:
    Array2d() : vector<T>(0) {}
    Array2d(size_t n1, size_t n2)
       : vector<T>(n1*n2), nrow(n1), ncol(n2) {}
    ~Array2d();
    
    T&       operator()(size_t i, size_t j)
	{ return vector<T>::operator[](offset(i,j)); }
    const T& operator()(size_t i, size_t j) const
	{ return vector<T>::operator[](offset(i,j)); }
    T* begin() { return &(vector<T>::front()); }
    T* end() { return (&(vector<T>::back())+1); }
    const T* begin() const { return &(vector<T>::front()); }
    const T* end() const { return (&(vector<T>::back())+1); }
    
    void resize(size_t n1, size_t n2, const T& val = T())
	{ vector<T>::resize(n1*n2, val); nrow=n1, ncol=n2; }

    size_t nCol() const { return ncol; }
    size_t nRow() const { return nrow; }

    // unary operators
    Array2d<T> operator-();
    
    // binary operators
    Array2d<T>& operator=(const T& val);
    Array2d<T>& operator+=(const Array2d<T>&);
    Array2d<T>& operator+=(const T& val);
    Array2d<T>& operator-=(const Array2d<T>&);
    Array2d<T>& operator-=(const T& val);
    Array2d<T>& operator*=(const T&);
    Array2d<T>& operator/=(const T&);
    
private:
    size_t size() const { return vector<T>::size(); }

    size_t offset(size_t i, size_t j) const;

    size_t nrow;
    size_t ncol;

};

template<class T> 
inline size_t Array2d<T>::offset(size_t i, size_t j) const
{
//    if ((i-1)*ncol+j-1 >=size()){
//	cerr << "Array2d OoR " << i << " " << j << " " << nrow << " " << ncol << '\n';
//    }
    return (i-1)*ncol+j-1;
}

template<class T>
Array2d<T>::~Array2d()
{
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest = val;
    return *this;
}

template<class T>
inline
Array2d<T> Array2d<T>::operator-()
{
    Array2d<T> negative(nrow,ncol);

    T* dest = negative.end();
    const T* src = end();
    while (dest > negative.begin()) *--dest = (-1)*(*--src);
    return negative;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator+=(const Array2d<T>& a)
{
    if (size() != a.size()) error("Array2d::operator+= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest += *--src;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator+=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest += val;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator-=(const Array2d<T>& a)
{
    if (size() != a.size()) error("Array2d::operator-= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest -= *--src;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator-=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest -= val;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator*=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest *= a;
    return *this;
}
    
template<class T>
inline
Array2d<T>& Array2d<T>::operator/=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest /= a;
    return *this;
}

// nonmember functions
template<class T>
inline
Array2d<T> operator+(const Array2d<T>& a, const Array2d<T>& b)
{
    Array2d<T> r=a;
    return r+=b;
}

template<class T>
inline
Array2d<T> operator-(const Array2d<T>& a, const Array2d<T>& b)
{
    Array2d<T> r=a;
    return r-=b;
}

template<class T>
inline
Array2d<T> operator*(const Array2d<T>& a, const T& val)
{
    Array2d<T> r=a;
    return r*=val;
}

template<class T>
inline
Array2d<T> operator*(const T& val, const Array2d<T>& a)
{
    return operator*(a,val);
}

template<class T>
inline
Array2d<T> operator/(const Array2d<T>& a, const T& val)
{
    Array2d<T> r=a;
    return r/=val;
}

template<class T>
ostream& operator<<(ostream& s, const Array2d<T>& a)
{
    for (int i=1; i<=a.nRow(); i++){
	for (int j=1; j<a.nCol(); j++){
	    s << a(i,j) << ", ";
	}
	s << a(i,a.nCol()) << '\n';
    }
    return s;
}

// multiplication
// note: this can be more optimized by directly accessing member data
//       if these are friend of Array class.
//       How to deal with offset() is the key?
//
// matrix * vector
template<class T>
Array1d<T> operator*(const Array2d<T>& A, const Array1d<T>& x)
{
    size_t m = A.nRow();
    size_t n = x.size();
    if (A.nCol() != n) error("matrix*vector::size mismatch");

    Array1d<T> b(m);
    for (size_t i=1; i<=m; i++){
	double val=0;
       for (size_t j=1; j<=n; j++){
	    val += A(i,j)*x(j);
	}
	b(i) = val;
    }

    return b;
}

// matrix * matrix
template<class T>
Array2d<T> operator*(const Array2d<T>& A, const Array2d<T>& B)
{
    size_t l = A.nRow();
    size_t m = A.nCol();
    size_t n = B.nCol();
    if (B.nRow() != m) error("matrix*matrix::size mismatch");

    Array2d<T> C(l,n);
    for (size_t i=1; i<=l; i++){
       for (size_t j=1; j<=n; j++){
	    double val=0;
           for (size_t k=1; k<=m; k++){
		val += A(i,k)*B(k,j);
	    }
	    C(i,j) = val;
	}
    }

    return C;
}

#endif /* _JK_NEW_ARRAY_H_ */
