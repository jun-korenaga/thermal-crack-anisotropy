/*
 * tensor.cc
 *
 * Summer 2024
 * Jun Korenaga
 */

#include <utility>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "emt/tensor.hpp"

using namespace std;

// for two-way Voigt indexing
int ij2I(int i, int j){
    static const int lookup[3][3] = {
        {1, 6, 5},
        {6, 2, 4},
        {5, 4, 3}
    };
    return lookup[i-1][j-1];
}

void I2ij(int I, int& i, int& j){
    static const int lookup_i[6] = {1, 2, 3, 2, 1, 1};
    static const int lookup_j[6] = {1, 2, 3, 3, 3, 2};
    i = lookup_i[I-1];
    j = lookup_j[I-1];
}

void Tensor4::initialize()
{
    A.resize(n,n);
    A = 0.0;
}

Tensor4::Tensor4()
{
    initialize();
}

Tensor4::Tensor4(Array1d<double>& c)
{
    initialize();

    switch(c.size()){
    case 2: // isotropic: c(1) = c11, c(2) = c44
    {
	double c11=c(1), c44=c(2), c12=c11-2.0*c44;
	A(1,1) = c11; A(1,2) = c12; A(1,3) = c12;
	A(2,1) = c12; A(2,2) = c11; A(2,3) = c12;
	A(3,1) = c12; A(3,2) = c12; A(3,3) = c11;
	A(4,4) = c44;
	A(5,5) = c44;
	A(6,6) = c44;
	break;
    }
    case 3: // cubic: c11, c12, c44
    {
	double c11=c(1), c12=c(2), c44=c(3);
	A(1,1) = c11; A(1,2) = c12; A(1,3) = c12;
	A(2,1) = c12; A(2,2) = c11; A(2,3) = c12;
	A(3,1) = c12; A(3,2) = c12; A(3,3) = c11;
	A(4,4) = c44;
	A(5,5) = c44;
	A(6,6) = c44;
	break;
    }
    case 5: // hexagonal (transversely isotropic): c11, c12, c13, c33, c44
    {
	double c11=c(1), c12=c(2), c13=c(3), c33=c(4), c44=c(5);
	double c66=0.5*(c11-c12);
	A(1,1) = c11; A(1,2) = c12; A(1,3) = c13;
	A(2,1) = c12; A(2,2) = c11; A(2,3) = c13;
	A(3,1) = c13; A(3,2) = c13; A(3,3) = c33;
	A(4,4) = c44;
	A(5,5) = c44;
	A(6,6) = c66;
	break;
    }
    case 9: // orthorhombic: c11, c12, c13, c22, c23, c33, c44, c55, c66
    {
	double c11=c(1), c12=c(2), c13=c(3), c22=c(4), c23=c(5);
	double c33=c(6), c44=c(7), c55=c(8), c66=c(9);
	A(1,1) = c11; A(1,2) = c12; A(1,3) = c13;
	A(2,1) = c12; A(2,2) = c22; A(2,3) = c23;
	A(3,1) = c13; A(3,2) = c23; A(3,3) = c33;
	A(4,4) = c44;
	A(5,5) = c55;
	A(6,6) = c66;
	break;
    }
    case 13: // monoclinic: c11, c12, c13, c15, c22, c23, c25, c33, c35,
             //             c44, c46, c55, c66
    {
	double c11=c(1), c12=c(2), c13=c(3), c15=c(4), c22=c(5);
	double c23=c(6), c25=c(7), c33=c(8), c35=c(9), c44=c(10);
	double c46=c(11), c55=c(12), c66=c(13);
	A(1,1) = c11; A(1,2) = c12; A(1,3) = c13; A(1,5)=c15;
	A(2,1) = c12; A(2,2) = c22; A(2,3) = c23; A(2,5)=c25;
	A(3,1) = c13; A(3,2) = c23; A(3,3) = c33; A(3,5)=c35;
	A(4,4) = c44; A(4,6) = c46;
	A(5,1) = c15; A(5,2) = c25; A(5,3) = c35; A(5,5) = c55;
	A(6,4) = c46; A(6,6) = c66;
	break;
    }
    default:
	cerr << "Tensor4: unknown input type detected.\n";
	exit(1);
	break;
    }
}

bool Tensor4::isMajorSymmetry() const
{
    bool isMajorSym=true;

    for (int i=1; i<=n; i++){
	for (int j=i+1; j<=n; j++){
	    if (A(i,j)!=A(j,i)){
		isMajorSym=false;
		break;
	    }
	}
    }

    return isMajorSym;
}

void Tensor4::setJ()
{
    // unit 4th-rank tensor
    // J_klij = (delta_ik * delta_jl + delta_il * delta_jk)/2
    A = 0.0;
    for (int k=1; k<=3; k++){
	for (int l=1; l<=3; l++){
	    for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++){
		    int delta_ik = (i==k) ? 1 : 0;
		    int delta_jl = (j==l) ? 1 : 0;
		    int delta_il = (i==l) ? 1 : 0;
		    int delta_jk = (j==k) ? 1 : 0;
		    A(ij2I(k,l),ij2I(i,j)) =
			0.5*(delta_ik * delta_jl + delta_il * delta_jk);
		}
	    }
	}
    }

    return;
}

void Tensor4::setTI(int i)
{
    // transversely isotropic basis tensors
    // kachanov18book eq. 1.4.13
    A = 0.0;
    
    switch(i){
    case 1:
	A(ij2I(1,1),ij2I(1,1)) = 1.0;
	A(ij2I(2,2),ij2I(2,2)) = 1.0;
	A(ij2I(1,1),ij2I(2,2)) = 1.0;
	A(ij2I(2,2),ij2I(1,1)) = 1.0;
	break;
    case 2:
	A(ij2I(1,2),ij2I(1,2)) = 0.5;
	A(ij2I(2,1),ij2I(2,1)) = 0.5;
	A(ij2I(1,2),ij2I(2,1)) = 0.5;
	A(ij2I(2,1),ij2I(1,2)) = 0.5;
	A(ij2I(1,1),ij2I(1,1)) = 0.5;
	A(ij2I(2,2),ij2I(2,2)) = 0.5;
	A(ij2I(1,1),ij2I(2,2)) = -0.5;
	A(ij2I(2,2),ij2I(1,1)) = -0.5;
	break;
    case 3:
	A(ij2I(1,1),ij2I(3,3)) = 1.0;
	A(ij2I(2,2),ij2I(3,3)) = 1.0;
	break;
    case 4:
	A(ij2I(3,3),ij2I(1,1)) = 1.0;
	A(ij2I(3,3),ij2I(2,2)) = 1.0;
	break;
    case 5:
	A(ij2I(1,3),ij2I(1,3)) = 0.25;
	A(ij2I(2,3),ij2I(2,3)) = 0.25;
	A(ij2I(1,3),ij2I(3,1)) = 0.25;
	A(ij2I(2,3),ij2I(3,2)) = 0.25;
	A(ij2I(3,1),ij2I(1,3)) = 0.25;
	A(ij2I(3,2),ij2I(2,3)) = 0.25;
	A(ij2I(3,1),ij2I(3,1)) = 0.25;
	A(ij2I(3,2),ij2I(3,2)) = 0.25;
	break;
    case 6:
	A(ij2I(3,3),ij2I(3,3)) = 1.0;
	break;
    default:
	cerr << "Tensor4:setTI() - unsupported index (" << i
	     << ") found.\n";
	exit(1);
    }
}


Tensor4 Tensor4::inverse() const
{
    Array2d<double> result(n,n), tmp(n,n);
    tmp = A; // create a copy of the original matrix

    // initialize result as identity matrix
    for (int i=1; i<=n; i++){
	result(i,i) = 1.0;
    }

    // perform Gauss-Jordan elimination
    for (int i=1; i<=n; i++){
	// find pivot
	int pivot=i;
	for (int j=i+1; j<=n; j++){
	    if (std::abs(tmp(j,i)) > std::abs(tmp(pivot,i))){
		pivot = j;
	    }
	}

	// swap rows if necessary
	if (pivot!=i){
	    for (int j=1; j<=n; j++){
		std::swap(tmp(i,j), tmp(pivot,j));
		std::swap(result(i,j), result(pivot,j));
	    }
	}

	// check for singularity
	if (tmp(i,i)==0.0){
	    cerr << "Tensor4::inverse() - singular and cannot be inverted.\n";
	    exit(1);
	}

	// scale the pivot row
	double scale = 1.0/tmp(i,i);
	for (int j=1; j<=n; j++){
	    tmp(i,j) *= scale;
	    result(i,j) *= scale;
	}

	// eliminate the pivot variable from other rows
	for (int j=1; j<=n; j++){
	    if (j!=i) {
		double factor = tmp(j, i);
		for (int k=1; k<=n; k++){
		    tmp(j,k) -= factor*tmp(i,k);
		    result(j,k) -= factor*result(i,k);
		}
	    }
	}
    }

    // adjustment for Voigt notation
    for (int i=1; i<=3; i++){
	for (int j=4; j<=6; j++){
	    result(i,j) = 0.5*result(i,j);
	    result(j,i) = result(i,j);
	}
    }
    for (int i=4; i<=6; i++){
	for (int j=4; j<=6; j++){
	    result(i,j) = 0.25*result(i,j);
	}
    }
    
    Tensor4 invA;
    invA.A = result;
    return invA;
}

Tensor4 Tensor4::operator-()
{
    Tensor4 negative;
    negative.A = -A;

    return negative;
}

Tensor4& Tensor4::operator=(double val)
{
    A = val;
    return *this;
}

Tensor4& Tensor4::operator+=(const Tensor4& B)
{
    A += B.A;
    return *this;
}

Tensor4& Tensor4::operator-=(const Tensor4& B)
{
    A -= B.A;
    return *this;
}

Tensor4& Tensor4::operator*=(double val)
{
    A *= val;
    return *this;
}

Tensor4& Tensor4::operator/=(double val)
{
    A /= val;
    return *this;
}

Tensor4 operator+(const Tensor4& A, const Tensor4& B)
{
    Tensor4 C;
    C.A = A.A + B.A;

    return C;
}

Tensor4 operator-(const Tensor4& A, const Tensor4& B)
{
    Tensor4 C;
    C.A = A.A - B.A;

    return C;
}

Tensor4 operator*(const Tensor4& A, const Tensor4& B)
{
    Tensor4 C;
    Array2d<double> tmpB(Tensor4::n,Tensor4::n);
    tmpB = B.A;
    for (int i=4; i<=Tensor4::n; i++){
	for (int j=1; j<=Tensor4::n; j++){
	    tmpB(i,j) *= 2.0;
	}
    }
    C.A = A.A * tmpB;

    return C;
}

bool operator==(const Tensor4& A, const Tensor4& B)
{
    bool is_equal=true;
    for (int i=1; i<=Tensor4::n; i++){
	for (int j=1; j<=Tensor4::n; j++){
	    if (A.A(i,j)!=B.A(i,j)){
		is_equal = false;
		break;
	    }
	}
    }

    return is_equal;
}

bool operator!=(const Tensor4& A, const Tensor4& B)
{
    return (A==B) ? false : true;
}

ostream& operator<<(ostream& os, const Tensor4& A)
{

    for (int i=1; i<=Tensor4::n; i++){
	for (int j=1; j<=Tensor4::n; j++){
	    double val = A.A(i,j);
	    if (std::abs(val)<1e-10) val = 0.0;

	    char buffer[50];
	    std::snprintf(buffer, sizeof(buffer), "%12.5g", val);

	    os << buffer << " ";
	}
	os << '\n';
    }
    return os;
}
