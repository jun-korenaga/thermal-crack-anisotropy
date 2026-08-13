/*
 * elasticity.h
 *
 * Summer 2024
 * Jun Korenaga
 */

#ifndef _JK_ELASTICITY_H_
#define _JK_ELASTICITY_H_

#include <array>
#include <fstream>
#include "emt/array.hpp"
#include "emt/ellipsoid.hpp"
#include "emt/tensor.hpp"

class HillTensor {
public:
    // default constructor
    HillTensor(const Tensor4& C, const Ellipsoid& e,
	       int ntheta, int nphi, double factor_w,
	       bool assume_iso=false,
	       int max_level=20,
	       bool dump_patches=false,
	       const std::string& dumpfn="");
    
    // constructor for simple recursive integration
    HillTensor(const Tensor4& C, const Ellipsoid& e,
	       int ntheta, int nphi, bool is_recursive=true,
	       int max_level=5, double recur_rel_tol=1e-5,
	       double recur_abs_tol=1e-8,
	       bool assume_iso=false);
    ~HillTensor() { if (p_dump!=nullptr) delete p_dump; }    

    HillTensor(const HillTensor&) = delete;
    HillTensor& operator=(const HillTensor&) = delete;

    void setTensor(Tensor4& _P){ _P = P; }
    int numIntegPoints(){ return n_int_points; }
    int numHitMaxLevel(){ return n_hit_max; }

private:
    static constexpr int sizeD=81;
    void set_U();
    void set_A();
    void set_B();
    void set_A_iso();
    void set_B_iso();
    double calc_f();
    void set_er();

    void set_D_flat();
    void set_D_recursive();
    void set_D_adaptive();
    std::array<double, sizeD> calc_dsum(double theta, double phi,
					double dtheta, double dphi);
    std::array<double, sizeD> calc_dsum_adaptive(double theta, double phi,
						 double dtheta, double dphi,
						 int level);
    std::array<double, sizeD> calc_dsum_recursive(double theta, double phi,
						  double dtheta, double dphi,
						  int level,
						  std::array<double, sizeD> old_dsum);
    Tensor4 C, P;
    Ellipsoid elp;
    Array2d<double> T;
    int ntheta, nphi, max_level;
    long int n_int_points;
    int n_hit_max;
    bool is_recursive;
    bool assume_iso;
    double f_min, factor_w, sin_theta_ref;
    double recur_rel_tol, recur_abs_tol;
    Array1d<int> diagD;

    Array1d<double> er, h_mj, A, dA_dtheta;
    Array2d<double> U, B, dB_dphi, dB_dtheta;
    Array2d< Array1d<double> > hmji;
    Array2d< Array2d<double> > D;
    Array2d< Array2d<int> > indexD;

    double sin_theta, cos_theta, sin_phi, cos_phi;
    double sin2_theta, cos2_theta, sin_cos_theta, sin4_theta, cos4_theta;
    double K, E, G, dK_dtheta, dE_dtheta, dG_dtheta;
    double c11, c12, c13, c33, c44, c66, L1, L2, L3;

    bool dump_patches;
    std::ofstream* p_dump;
};

Tensor4 eshelby_isotropic_exact(const Tensor4& C, double a1, double a3);

Tensor4 eshelby_transverse_exact(const Tensor4& C, double a1, double a3);

#endif /* _JK_ELASTICITY_H_ */
