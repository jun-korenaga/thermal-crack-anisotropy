/*
 * ellipsoid.cc
 *
 * Summer 2024
 * Jun Korenaga
 */

#include <cmath>
#include "emt/array.hpp"
#include "emt/ellipsoid.hpp"

using namespace std;

Ellipsoid::Ellipsoid(double _a1, double _a2, double _a3,
		     double _alpha, double _beta, double _gamma)
    : a1(_a1), a2(_a2), a3(_a3),
      alpha(_alpha), beta(_beta), gamma(_gamma)
{
    T.resize(3,3);
    T = 0.0;
    T(1,1) = 1.0/(a1*a1);
    T(2,2) = 1.0/(a2*a2);
    T(3,3) = 1.0/(a3*a3);

    if (alpha!=0.0 || beta!=0.0 || gamma!=0.0){
	// needs rotation

	double deg2rad = M_PI/180.0;
	alpha *= deg2rad;
	beta *= deg2rad;
	gamma *= deg2rad;
	
	Array2d<double> R(3,3), Rt(3,3);
	double cos_alpha=cos(alpha), sin_alpha=sin(alpha);
	double cos_beta=cos(beta), sin_beta=sin(beta);
	double cos_gamma=cos(gamma), sin_gamma=sin(gamma);
	R(1,1) = cos_alpha*cos_beta;
	R(1,2) = cos_alpha*sin_beta*sin_gamma - sin_alpha*cos_gamma;
	R(1,3) = cos_alpha*sin_beta*cos_gamma + sin_alpha*sin_gamma;
	R(2,1) = sin_alpha*cos_beta;
	R(2,2) = sin_alpha*sin_beta*sin_gamma + cos_alpha*cos_gamma;
	R(2,3) = sin_alpha*sin_beta*cos_gamma - cos_alpha*sin_gamma;
	R(3,1) = -sin_beta;
	R(3,2) = cos_beta*sin_gamma;
	R(3,3) = cos_beta*cos_gamma;

	for (int i=1; i<=3; i++){
	    for (int j=1; j<=3; j++){
		Rt(i,j) = R(j,i); // transpose of R
	    }
	}
	T = R*T;
	T = T*Rt;
    }
}

void Ellipsoid::setTensor(Array2d<double>& _T) const
{
    if (_T.nRow()==3 && _T.nCol()==3){
	for (int i=1; i<=3; i++){
	    for (int j=1; j<=3; j++){
		_T(i,j) = T(i,j);
	    }
	}
    }else{
	cerr << "Ellipsoid::setTensor - incorrect input array dimension\n";
	exit(1);
    }
}
