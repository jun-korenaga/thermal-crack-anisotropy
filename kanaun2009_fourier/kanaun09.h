/*
 * kanaun09.h
 *
 * Fourier-domain solution for an elliptical crack under a constant
 * external stress field, based on Kanaun and Levin (2009), Eqs. (23)-(32).
 *
 * Summer 2026
 * Jun Korenaga
 */

#ifndef _JK_KANAUN09_H_
#define _JK_KANAUN09_H_

#include "emt/tensor.hpp"

class Kanaun09 {
public:
    Kanaun09(const Tensor4& C, double a1, double a2,
	      double alpha=0.0, double beta=0.0, double gamma=0.0,
	      int nphi=360, int nq=64);

    void setT0(Array2d<double>& _T0) const { _T0 = T0; }
    void setM(Array2d<double>& _M) const { _M = M; }
    void setCrackCompliance(Tensor4& _H) const { _H = H; }
    long int numIntegPoints() const { return n_int_points; }

private:
    void set_R();
    void set_C_local();
    void set_quadrature();
    void set_T0();
    void set_M();
    void set_H();
    void calc_S(double q1, double q2, double q3,
		Array2d<double>& S) const;
    void invert3(const Array2d<double>& A, Array2d<double>& invA) const;

    Tensor4 C, C_local, H;
    Array2d<double> R, T0, M;
    Array1d<double> xq, wq;
    double a1, a2, alpha, beta, gamma;
    int nphi, nq;
    long int n_int_points;
};

void kanaun09_isotropic_T0(const Tensor4& C, double a,
			   Array2d<double>& T0);
void kanaun09_transverse_T0(const Tensor4& C, double a,
			    Array2d<double>& T0);

#endif /* _JK_KANAUN09_H_ */
