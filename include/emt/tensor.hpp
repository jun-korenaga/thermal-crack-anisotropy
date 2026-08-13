/*
 * tensor.h
 *
 * Summer 2024
 * Jun Korenaga
 */

#ifndef _JK_TENSOR_H_
#define _JK_TENSOR_H_

#include <iostream>
#include "emt/array.hpp"

// for two-way Voigt indexing
int ij2I(int i, int j);
void I2ij(int I, int& i, int& j);

// 4th-rank tensor class
class Tensor4 {
public:
    Tensor4();

    // the length of the array c can be 2 (isotropic, c11, c44),
    // 3 (cubic, c11, c12, c44)
    // 5 (hexagonal; transversely isotropic, c11, c12, c13, c33, c44),
    // 9 (orthorhombic, c11, c12, c13, c22, c23, c33, c44, c55, c66),
    // and 13 (monoclinic, c11, c12, c13, c15, c22, c23, c25, c33, c35,
    // c44, c46, c55, c66)
    Tensor4(Array1d<double>& c); 

    bool isMajorSymmetry() const; 
    double val(int i, int j, int k, int l) const
	{ return A(ij2I(i,j),ij2I(k,l)); }
    double val(int I, int J) const
	{ return A(I,J); }
    void set(int i, int j, int k, int l, double val)
	{ A(ij2I(i,j),ij2I(k,l))=val; return; }
    void set(int I, int J, double val)
	{ A(I,J)=val; return; }
    void set(Array2d<double>& _A){ A = _A; }

    void setJ(); // J_ijkl = (delta_ik * delta_jl + delta_il * delta_jk)/2
    void setTI(int i); // i-th transversely isotropic basis tensor
    Tensor4 inverse() const;

    // unary operators
    Tensor4 operator-();
    
    // binary operators
    Tensor4& operator=(double val);
    Tensor4& operator+=(const Tensor4& B);
    Tensor4& operator-=(const Tensor4& B);
    Tensor4& operator*=(double val);
    Tensor4& operator/=(double val);

    // helper functions
    friend Tensor4 operator+(const Tensor4& A, const Tensor4& B);
    friend Tensor4 operator-(const Tensor4& A, const Tensor4& B);
    friend Tensor4 operator*(const Tensor4& A, const Tensor4& B);
    friend bool operator==(const Tensor4& A, const Tensor4& B);
    friend bool operator!=(const Tensor4& A, const Tensor4& B);
    friend std::ostream& operator<<(std::ostream&, const Tensor4&);
    
private:
    static constexpr int n=6;
    void initialize();

    // for Voigt indexing
    Array2d<double> A;
};

// Tensor4 + Tensor4
Tensor4 operator+(const Tensor4& A, const Tensor4& B);

// Tensor4 - Tensor4
Tensor4 operator-(const Tensor4& A, const Tensor4& B);

// Tensor4 * Tensor4
Tensor4 operator*(const Tensor4& A, const Tensor4& B);

// comparison operators
bool operator==(const Tensor4& A, const Tensor4& B);
bool operator!=(const Tensor4& A, const Tensor4& B);

// output
std::ostream& operator<<(std::ostream&, const Tensor4&);

#endif /* _JK_TENSOR_H_ */
