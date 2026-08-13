/*
 * hill_tensor.c
 *
 * Summer 2024
 * Jun Korenaga
 */

#include "emt/hill_tensor.hpp"

using namespace std;

HillTensor::HillTensor(const Tensor4& _C, const Ellipsoid& _e,
		       int _ntheta, int _nphi,
		       double _w, bool _iso, int _mlevel,
		       bool dump_patches, const std::string& dumpfn)
    : C(_C), elp(_e), ntheta(_ntheta), nphi(_nphi),
      max_level(_mlevel), n_int_points(0), n_hit_max(0),
      is_recursive(false), assume_iso(_iso), f_min(0.0), factor_w(_w),
      sin_theta_ref(0.0), recur_rel_tol(0.0), recur_abs_tol(0.0),
      dump_patches(dump_patches), p_dump(nullptr)
{
    if (dump_patches && !dumpfn.empty()) {
	p_dump = new ofstream(dumpfn);
	if (!p_dump->is_open()) {
	    std::cerr << "Failed to open file: " << dumpfn << endl;
	    delete p_dump;
	    p_dump = nullptr;
	}
    }
    
    // set the ellipsoid tensor
    T.resize(3,3);
    elp.setTensor(T);

    // prepare radial unit vector
    er.resize(3);

    // set stiffness-related constants
    c11=C.val(1,1); c12=C.val(1,2); c13=C.val(1,3);
    c33=C.val(3,3); c44=C.val(4,4); c66=C.val(6,6);
    L1 = c44/c66;
    L2 = (c11*c33-c13*c13-2*c13*c44)/(c11*c44);
    L3 = c33/c11;

    // set up other temporary arrays
    D.resize(3,3);
    indexD.resize(3,3);
    int iD=0;
    for (int i=1; i<=3; i++){
	for (int j=1; j<=3; j++){
	    D(i,j).resize(3,3);
	    D(i,j) = 0.0;

	    indexD(i,j).resize(3,3);
	    for (int k=1; k<=3; k++){
		for (int l=1; l<=3; l++){
		    indexD(i,j)(k,l) = iD++;
		}
	    }
	}
    }
    h_mj.resize(3);
    U.resize(3,3);
    A.resize(3); dA_dtheta.resize(3);
    B.resize(3,3); dB_dphi.resize(3,3); dB_dtheta.resize(3,3);

    set_D_adaptive();
    for (int I=1; I<=6; I++){
	int m, p;
	I2ij(I,m,p);
	for (int J=1; J<=6; J++){
	    int i, j;
	    I2ij(J,i,j);

	    double val = 0.25*(D(m,j)(i,p) + D(p,j)(i,m)
			       + D(m,i)(j,p) + D(p,i)(j,m));
	    val *= -1.0/(4*M_PI); 
	    P.set(I,J,val);
	}
    }
}

HillTensor::HillTensor(const Tensor4& _C, const Ellipsoid& _e,
		       int _ntheta, int _nphi,
		       bool _recur, int _max_level,
		       double _rel_tol, double _abs_tol,
		       bool _iso)
    : C(_C), elp(_e), ntheta(_ntheta), nphi(_nphi),
      max_level(_max_level), n_int_points(0), n_hit_max(0),
      is_recursive(_recur), assume_iso(_iso), f_min(0.0), factor_w(0.0),
      sin_theta_ref(0.0), recur_rel_tol(_rel_tol), recur_abs_tol(_abs_tol),
      dump_patches(false), p_dump(nullptr)
{
    // set the ellipsoid tensor
    T.resize(3,3);
    elp.setTensor(T);

    // prepare radial unit vector
    er.resize(3);

    // set stiffness-related constants
    c11=C.val(1,1); c12=C.val(1,2); c13=C.val(1,3);
    c33=C.val(3,3); c44=C.val(4,4); c66=C.val(6,6);
    L1 = c44/c66;
    L2 = (c11*c33-c13*c13-2*c13*c44)/(c11*c44);
    L3 = c33/c11;

    // set up other temporary arrays
    D.resize(3,3);
    indexD.resize(3,3);
    diagD.resize(6);
    int iD=0;
    int idiag=1;
    for (int i=1; i<=3; i++){
	for (int j=1; j<=3; j++){
	    D(i,j).resize(3,3);
	    D(i,j) = 0.0;

	    indexD(i,j).resize(3,3);
	    for (int k=1; k<=3; k++){
		for (int l=1; l<=3; l++){
		    if (i==j && j==k && k==l){
			diagD(idiag++) = iD;
		    }
		    indexD(i,j)(k,l) = iD++;
		}
	    }
	}
    }
    h_mj.resize(3);
    U.resize(3,3);
    A.resize(3); dA_dtheta.resize(3);
    B.resize(3,3); dB_dphi.resize(3,3); dB_dtheta.resize(3,3);

    if (is_recursive){
	set_D_recursive();
    }else{
	set_D_flat();
    }
    for (int I=1; I<=6; I++){
	int m, p;
	I2ij(I,m,p);
	for (int J=1; J<=6; J++){
	    int i, j;
	    I2ij(J,i,j);

	    double val = 0.25*(D(m,j)(i,p) + D(p,j)(i,m)
			       + D(m,i)(j,p) + D(p,i)(j,m));
	    val *= -1.0/(4*M_PI); 
	    P.set(I,J,val);
	}
    }
}

void HillTensor::set_D_adaptive()
{
    std::array<double, sizeD> dsum;
    dsum.fill(0.0);
    
    // set up initial integration points (centroids) on a sphere
    // theta: 0 to PI (note: the centroid approach naturally avoids 0 and PI)
    // phi: 0 to 2*PI
    Array1d<double> theta(ntheta), phi(nphi);
    double dtheta = M_PI/ntheta; 
    double dphi = 2*M_PI/nphi;
    for (int ii=1; ii<=ntheta; ii++){
	theta(ii) = dtheta*(ii-0.5);
    }
    for (int ii=1; ii<=nphi; ii++){
	phi(ii) = dphi*(ii-0.5);
    }

    // adaoptive integration
    n_int_points = 0;
    n_hit_max = 0;
    double a1, a2, a3;
    elp.setSemiAxes(a1,a2,a3);
    double a_min = std::min({a1,a2,a3});
    f_min = a_min*a_min;
    sin_theta_ref = sin(M_PI*0.25);

    for (int ii=1; ii<=ntheta; ii++){
	for (int jj=1; jj<=nphi; jj++){
	    dsum = calc_dsum_adaptive(theta(ii),phi(jj),dtheta,dphi,0);

	    for (int m=1; m<=3; m++){
		for (int j=1; j<=3; j++){
		    for (int i=1; i<=3; i++){
			for (int p=1; p<=3; p++){
			    D(m,j)(i,p) += dsum[indexD(m,j)(i,p)];
			} // for p
		    } // for i
		} // for j
	    } // for m
	    
	} // for jj
    } // for ii
}

void HillTensor::set_D_flat()
{
    n_int_points = ntheta*nphi;
    n_hit_max = 0;
    
    // set up integration points (centroids) on a sphere
    // theta: 0 to PI (note: the centroid approach naturally avoids 0 and PI)
    // phi: 0 to 2*PI
    Array1d<double> theta(ntheta), phi(nphi);
    double dtheta = M_PI/ntheta; 
    double dphi = 2*M_PI/nphi;
    for (int ii=1; ii<=ntheta; ii++){
	theta(ii) = dtheta*(ii-0.5);
    }
    for (int ii=1; ii<=nphi; ii++){
	phi(ii) = dphi*(ii-0.5);
    }

    // centroid integration 
    for (int ii=1; ii<=ntheta; ii++){
	sin_theta = sin(theta(ii));
	cos_theta = cos(theta(ii));
	double dS = sin_theta*dtheta*dphi;
	if (assume_iso){
	    set_A_iso();
	}else{
	    set_A();
	}
	    
	for (int jj=1; jj<=nphi; jj++){
	    sin_phi = sin(phi(jj));
	    cos_phi = cos(phi(jj));

	    set_U();
	    if (assume_iso){
		set_B_iso();
	    }else{
		set_B();
	    }
	    set_er();
	    double f = calc_f();

	    for (int m=1; m<=3; m++){
		for (int j=1; j<=3; j++){
		    int k = (m==3 || j==3) ? 2 : 1;
		    h_mj(1) = -A(k)*B(m,j);
		    h_mj(2) = dA_dtheta(k)*B(m,j) + A(k)*dB_dtheta(m,j);
		    h_mj(3) = A(k)*dB_dphi(m,j);

		    for (int i=1; i<=3; i++){
			double Emji=0.0;
			for (int q=1; q<=3; q++){
			    Emji += U(i,q)*h_mj(q); 
			}

			for (int p=1; p<=3; p++){
			    double sum=0.0;
			    for (int s=1; s<=3; s++){
				sum += T(p,s)*er(s);
			    }
			    D(m,j)(i,p) += sum*f*Emji*dS;
			} // for p
		    } // for i
		} // for j
	    } // for m
	} // for jj
    } // for ii
}

void HillTensor::set_D_recursive()
{
    std::array<double, sizeD> dsum;
    dsum.fill(0.0);
    
    // recursive integration
    n_int_points = 0;
    n_hit_max = 0;
    
    // set up initial integration points (centroids) on a sphere
    // theta: 0 to PI (note: the centroid approach naturally avoids 0 and PI)
    // phi: 0 to 2*PI
    Array1d<double> theta(ntheta), phi(nphi);
    double dtheta = M_PI/ntheta; 
    double dphi = 2*M_PI/nphi;
    for (int ii=1; ii<=ntheta; ii++){
	theta(ii) = dtheta*(ii-0.5);
    }
    for (int ii=1; ii<=nphi; ii++){
	phi(ii) = dphi*(ii-0.5);
    }

    // recursive integration 
    for (int ii=1; ii<=ntheta; ii++){
	for (int jj=1; jj<=nphi; jj++){
	    dsum = calc_dsum_recursive(theta(ii),phi(jj),dtheta,dphi,0,dsum);

	    for (int m=1; m<=3; m++){
		for (int j=1; j<=3; j++){
		    for (int i=1; i<=3; i++){
			for (int p=1; p<=3; p++){
			    D(m,j)(i,p) += dsum[indexD(m,j)(i,p)];
			} // for p
		    } // for i
		} // for j
	    } // for m
	    
	} // for jj
    } // for ii
}

std::array<double, HillTensor::sizeD> HillTensor::calc_dsum(double theta, double phi,
							    double dtheta, double dphi)
{
    if (dump_patches && p_dump != nullptr){
    	*p_dump << theta << " " << phi << " " << dtheta << " " << dphi << endl;
    }
    
    std::array<double, sizeD> dsum;

    n_int_points++;

    // calculate with a given patch
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    sin_phi = sin(phi);
    cos_phi = cos(phi);

    // trapezoidal rule here
    double dS = 0.5*(sin(theta-0.5*dtheta)+sin(theta+0.5*dtheta))*dtheta*dphi;

    // note: I hope the use of global variables (cos2_theta, A, B, etc.) is OK...
    if (assume_iso){
	set_A_iso();
	set_B_iso();
    }else{
	set_A();
	set_B();
    }
    set_U();
    set_er();
    double f = calc_f();

    for (int m=1; m<=3; m++){
	for (int j=1; j<=3; j++){
	    int k = (m==3 || j==3) ? 2 : 1;
	    h_mj(1) = -A(k)*B(m,j);
	    h_mj(2) = dA_dtheta(k)*B(m,j) + A(k)*dB_dtheta(m,j);
	    h_mj(3) = A(k)*dB_dphi(m,j);
	    
	    for (int i=1; i<=3; i++){
		double Emji=0.0;
		for (int q=1; q<=3; q++){
		    Emji += U(i,q)*h_mj(q); 
		}

		for (int p=1; p<=3; p++){
		    double sum=0.0;
		    for (int s=1; s<=3; s++){
			sum += T(p,s)*er(s);
		    }
		    dsum[indexD(m,j)(i,p)] = sum*f*Emji*dS;
		} // for p
	    } // for i
	} // for j
    } // for m

    return dsum;
}

std::array<double, HillTensor::sizeD>
HillTensor::calc_dsum_adaptive(double theta, double phi,
			       double dtheta, double dphi,
			       int level)
{
    std::array<double, sizeD> dsum, dsum1, dsum2, dsum3, dsum4;
    
    // check factor f
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    set_er();
    double f = calc_f();
    int target_level = std::round(1.0-sin_theta_ref/sin_theta+factor_w*log(f/f_min));
    if (target_level<0) target_level = 0;
    if (target_level>max_level){
	target_level = max_level;
	n_hit_max++;
    }

    if (level>max_level+10){
	cerr << "theta=" << theta << ", phi= " << phi << ", level=" << level << '\n';
	cerr << "target_level=" << target_level << '\n';
	exit(1);
    }

    if (level>=target_level){
	// note: level>target_level can happen during mesh refinement
	dsum = calc_dsum(theta,phi,dtheta,dphi);
    }else{
	// make finer patches
	double dtheta05 = dtheta*0.5, dphi05 = dphi*0.5;
	double dtheta025 = dtheta*0.25, dphi025 = dphi*0.25;
	int level1 = level+1;

	dsum1 = calc_dsum_adaptive(theta-dtheta025,phi-dphi025,dtheta05,dphi05,level1);
	dsum2 = calc_dsum_adaptive(theta+dtheta025,phi-dphi025,dtheta05,dphi05,level1);
	dsum3 = calc_dsum_adaptive(theta-dtheta025,phi+dphi025,dtheta05,dphi05,level1);
	dsum4 = calc_dsum_adaptive(theta+dtheta025,phi+dphi025,dtheta05,dphi05,level1);
	for (int i=0; i<sizeD; i++){
	    dsum[i] = dsum1[i] + dsum2[i] + dsum3[i] + dsum4[i];
	}
    }
    return dsum;
}

std::array<double, HillTensor::sizeD>
HillTensor::calc_dsum_recursive(double theta, double phi,
				double dtheta, double dphi,
				int level,
				std::array<double, sizeD> old_dsum)
{
    std::array<double, sizeD> new_dsum, dsum1, dsum2, dsum3, dsum4;

    // make finer patches
    double dtheta05 = dtheta*0.5, dphi05 = dphi*0.5;
    double dtheta025 = dtheta*0.25, dphi025 = dphi*0.25;
    int level1 = level+1;

    dsum1 = calc_dsum(theta-dtheta025,phi-dphi025,dtheta05,dphi05);
    dsum2 = calc_dsum(theta+dtheta025,phi-dphi025,dtheta05,dphi05);
    dsum3 = calc_dsum(theta-dtheta025,phi+dphi025,dtheta05,dphi05);
    dsum4 = calc_dsum(theta+dtheta025,phi+dphi025,dtheta05,dphi05);
    for (int i=0; i<sizeD; i++){
	new_dsum[i] = dsum1[i] + dsum2[i] + dsum3[i] + dsum4[i];
    }
    
    // check convergence
    double dsum_sq=0.0, diff_sq=0.0;
    for (int i=1; i<=6; i++){
	int j = diagD(i);
    	double diff = new_dsum[j]-old_dsum[j];
	dsum_sq += old_dsum[j]*old_dsum[j];
	diff_sq += diff*diff;
    }

    double rel_error = sqrt(diff_sq/(dsum_sq+1e-15));
    double abs_error = sqrt(diff_sq);
    if (rel_error<recur_rel_tol || abs_error<recur_abs_tol || level1>max_level){
	if (rel_error>=recur_rel_tol && abs_error>=recur_abs_tol
	    && level1>max_level) n_hit_max++;
	return new_dsum;
    }else{
	dsum1 = calc_dsum_recursive(theta-dtheta025,phi-dphi025,
				    dtheta05,dphi05,level1,dsum1);
	dsum2 = calc_dsum_recursive(theta+dtheta025,phi-dphi025,
				    dtheta05,dphi05,level1,dsum2);
	dsum3 = calc_dsum_recursive(theta-dtheta025,phi+dphi025,
				    dtheta05,dphi05,level1,dsum3);
	dsum4 = calc_dsum_recursive(theta+dtheta025,phi+dphi025,
				    dtheta05,dphi05,level1,dsum4);
	for (int i=0; i<sizeD; i++){
	    new_dsum[i] = dsum1[i] + dsum2[i] + dsum3[i] + dsum4[i];
	}
	return new_dsum;
    }
}

void HillTensor::set_er()
{
    // radial unit vector
    er(1) = sin_theta*cos_phi;
    er(2) = sin_theta*sin_phi;
    er(3) = cos_theta;
}

double HillTensor::calc_f()
{
    double sum=0.0;
    for (int i=1; i<=3; i++){
	for (int j=1; j<=3; j++){
	    sum += er(i)*T(i,j)*er(j);
	}
    }

    return 1.0/sum;
}

void HillTensor::set_U()
{
    U(1,1) = sin_theta*cos_phi;
    U(1,2) = cos_theta*cos_phi;
    U(1,3) = -sin_phi/sin_theta;
    U(2,1) = sin_theta*sin_phi;
    U(2,2) = cos_theta*sin_phi;
    U(2,3) = cos_phi/sin_theta;	
    U(3,1) = cos_theta;
    U(3,2) = -sin_theta;
    U(3,3) = 0.0;
}

void HillTensor::set_A_iso()
{
    A(1) = 1.0/(2.0*c11*c66);
    A(2) = A(1);
    dA_dtheta(1) = 0.0;
    dA_dtheta(2) = 0.0;

    cos2_theta = cos_theta*cos_theta;
    sin2_theta = sin_theta*sin_theta;
    sin_cos_theta = sin_theta*cos_theta;
}

void HillTensor::set_A()
{
    cos2_theta = cos_theta*cos_theta;
    sin2_theta = sin_theta*sin_theta;
    sin_cos_theta = sin_theta*cos_theta;
    cos4_theta = cos2_theta*cos2_theta;
    sin4_theta = sin2_theta*sin2_theta;

    K = sqrt(cos2_theta+L1*sin2_theta);
    E = sqrt(cos4_theta+L2*sin2_theta*cos2_theta+L3*sin4_theta);
    G = sqrt(2.0*cos2_theta+L2*sin2_theta+2*E);

    dK_dtheta = sin_theta*cos_theta*(L1-1)/K;
    dE_dtheta = (-2.0*sin_cos_theta*cos2_theta
		 + L2*sin_theta*cos_theta*(cos2_theta-sin2_theta)
		 + 2.0*L3*sin_cos_theta*sin2_theta)/E;
    dG_dtheta = ((L2-2.0)*sin_cos_theta+dE_dtheta)/G;

    A(1) = 1.0/(c11*c44*K*E*G*sin2_theta);
    A(2) = 1.0/(c11*c44*E*G);
    dA_dtheta(1) = -A(1)*(dK_dtheta/K + dE_dtheta/E
			 + dG_dtheta/G + 2*cos_theta/sin_theta);
    dA_dtheta(2) = -A(2)*(dE_dtheta/E + dG_dtheta/G);
}

void HillTensor::set_B_iso()
{
    double cos2_phi = cos_phi*cos_phi;
    double sin2_phi = sin_phi*sin_phi;
    double sin_cos_phi = sin_phi*cos_phi;

    B(1,1) = c11+c66+(c11-c66)*sin2_theta*cos2_phi;
    B(1,2) = (c12+c66)*sin2_theta*sin_cos_phi;
    B(1,3) = (c12+c66)*sin_cos_theta*cos_phi;
    B(2,2) = c11+c66+(c11-c66)*sin2_theta*sin2_phi;
    B(2,3) = (c12+c66)*sin_cos_theta*sin_phi;
    B(3,3) = c11+c11*cos2_theta+c66*sin2_theta;

    dB_dphi(1,1) = (c11-c66)*sin2_theta*(-2.0*sin_cos_phi);
    dB_dphi(1,2) = (c12+c66)*sin2_theta*(cos2_phi-sin2_phi);
    dB_dphi(1,3) = (c12+c66)*sin_cos_theta*(-sin_phi);
    dB_dphi(2,2) = (c11-c66)*sin2_theta*(2.0*sin_cos_phi);
    dB_dphi(2,3) = (c12+c66)*sin_cos_theta*cos_phi;
    dB_dphi(3,3) = 0.0;

    dB_dtheta(1,1) = (c11-c66)*(2.0*sin_cos_theta)*cos2_phi;
    dB_dtheta(1,2) = (c12+c66)*(2.0*sin_cos_theta)*sin_cos_phi;
    dB_dtheta(1,3) = (c12+c66)*(cos2_theta-sin2_theta)*cos_phi;
    dB_dtheta(2,2) = (c11-c66)*(2.0*sin_cos_theta)*sin2_phi;
    dB_dtheta(2,3) = (c12+c66)*(cos2_theta-sin2_theta)*sin_phi;
    dB_dtheta(3,3) = c11*(-2.0*sin_cos_theta)+c66*(2.0*sin_cos_theta);

    for (int i=2; i<=3; i++){
	for (int j=1; j<i; j++){
	    B(i,j) = B(j,i);
	    dB_dphi(i,j) = dB_dphi(j,i);
	    dB_dtheta(i,j) = dB_dtheta(j,i);
	}
    }
}

void HillTensor::set_B()
{
    double cos2_phi = cos_phi*cos_phi;
    double sin2_phi = sin_phi*sin_phi;
    double sin_cos_phi = sin_phi*cos_phi;

    double xi1 = K*(c11*E+c11*cos2_theta+c44*sin2_theta)-c11*G*cos2_theta;
    double xi2 = c11*G*(cos4_theta+K*K*E-K*G*cos2_theta);
    double dxi1_dtheta = dK_dtheta*(c11*E+c11*cos2_theta+c44*sin2_theta)
	+ K*(c11*dE_dtheta-2.0*c11*sin_cos_theta+2.0*c44*sin_cos_theta)
	- c11*dG_dtheta*cos2_theta + 2.0*c11*G*sin_cos_theta;
    double dxi2_dtheta = c11*dG_dtheta*(cos4_theta+K*K*E-K*G*cos2_theta)
	+ c11*G*(-4.0*sin_cos_theta*cos2_theta + 2.0*K*E*dK_dtheta
		 + K*K*dE_dtheta - G*dK_dtheta*cos2_theta
		 - K*dG_dtheta*cos2_theta + 2.0*K*G*sin_cos_theta);

    // B
    B(1,1) = (cos2_theta*cos2_phi + E*sin2_phi)*xi1 + xi2*cos2_phi;
    B(1,2) = ((cos2_theta-E)*xi1+xi2)*sin_cos_phi;
    B(1,3) = (c13+c44)*sin_cos_theta*cos_phi;
    B(2,2) = (cos2_theta*sin2_phi + E*cos2_phi)*xi1 + xi2*sin2_phi;
    B(2,3) = (c13+c44)*sin_cos_theta*sin_phi;
    B(3,3) = c11*E+c11*cos2_theta+c44*sin2_theta;

    // dB_dphi
    dB_dphi(1,1) = 2.0*sin_cos_phi*((E-cos2_theta)*xi1-xi2);
    dB_dphi(1,2) = ((cos2_theta-E)*xi1+xi2)*(cos2_phi-sin2_phi);
    dB_dphi(1,3) = -(c13+c44)*sin_cos_theta*sin_phi;
    dB_dphi(2,2) = 2.0*sin_cos_phi*((cos2_theta-E)*xi1+xi2);
    dB_dphi(2,3) = (c13+c44)*sin_cos_theta*cos_phi;
    dB_dphi(3,3) = 0.0;

    // dB_dtheta
    dB_dtheta(1,1) = (-2.0*sin_cos_theta*cos2_phi + dE_dtheta*sin2_phi)*xi1
	+ (cos2_theta*cos2_phi + E*sin2_phi)*dxi1_dtheta
	+ cos2_phi*dxi2_dtheta;
    dB_dtheta(1,2) = ((-2.0*sin_cos_theta-dE_dtheta)*xi1
		      + (cos2_theta-E)*dxi1_dtheta + dxi2_dtheta)*sin_cos_phi;
    dB_dtheta(1,3) = (c13+c44)*(cos2_theta-sin2_theta)*cos_phi;
    dB_dtheta(2,2) = (-2.0*sin_cos_theta*sin2_phi + dE_dtheta*cos2_phi)*xi1
	+ (cos2_theta*sin2_phi + E*cos2_phi)*dxi1_dtheta
	+ sin2_phi*dxi2_dtheta;
    dB_dtheta(2,3) = (c13+c44)*(cos2_theta-sin2_theta)*sin_phi;
    dB_dtheta(3,3) = c11*dE_dtheta + 2.0*sin_cos_theta*(c44-c11);

    for (int i=2; i<=3; i++){
	for (int j=1; j<i; j++){
	    B(i,j) = B(j,i);
	    dB_dphi(i,j) = dB_dphi(j,i);
	    dB_dtheta(i,j) = dB_dtheta(j,i);
	}
    }
}
