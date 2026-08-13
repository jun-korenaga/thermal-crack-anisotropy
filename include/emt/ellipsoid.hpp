/*
 * ellipsoid.h
 *
 * Summer 2024
 * Jun Korenaga
 */

#ifndef _JK_ELLIPSOID_H_
#define _JK_ELLIPSOID_H_

#include "emt/array.hpp"

class Ellipsoid {
public:
    Ellipsoid(double _a1, double _a2, double _a3,
	      double _alpha=0.0, double _beta=0.0, double _gamma=0.0);
    // a1, a2, a3: semi-axes 
    // alpha (yaw), beta (pitch), gamma (roll):
    // intrinsic Euler angles about z-y'-x'' (in degree)

    void setTensor(Array2d<double>& _T) const;
    void setSemiAxes(double& _a1, double& _a2, double& _a3)
	{ _a1 = a1; _a2 = a2; _a3 = a3; }
private:
    double a1, a2, a3;
    double alpha, beta, gamma;
    Array2d<double> T;
};


#endif /* _JK_ELLIPSOID_H_ */
