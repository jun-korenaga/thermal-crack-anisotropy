/*
 * benchmark_kanaun09.cc
 *
 * usage: benchmark_kanaun09 [ -C<c11/c12/c13/c33/c44> ]
 *                           [ -P<nphi> -Q<nq> -N<ntheta> -W<factor_w> ]
 *                           [ -F ]
 *
 * -F: skip the finite-aspect-ratio comparison
 *
 * Summer 2026
 * Jun Korenaga
 */

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "emt/effective_medium.hpp"
#include "emt/ellipsoid.hpp"
#include "emt/hill_tensor.hpp"
#include "kanaun09.h"

using namespace std;

double relative_rms(const Tensor4& A, const Tensor4& B)
{
    double diff_sq=0.0, ref_sq=0.0;
    for (int I=1; I<=6; I++){
	for (int J=1; J<=6; J++){
	    double diff=A.val(I,J)-B.val(I,J);
	    diff_sq += diff*diff;
	    ref_sq += B.val(I,J)*B.val(I,J);
	}
    }
    return sqrt(diff_sq/ref_sq);
}

double max_abs_diff(const Tensor4& A, const Tensor4& B)
{
    double max_diff=0.0;
    for (int I=1; I<=6; I++){
	for (int J=1; J<=6; J++){
	    double diff=abs(A.val(I,J)-B.val(I,J));
	    if (diff>max_diff) max_diff=diff;
	}
    }
    return max_diff;
}

void print_row(const string& section, const string& test,
	       double dip, double aspect, const string& component,
	       double value, double reference, double error,
	       int nphi, int nq, int ntheta,
	       long int nk, long int na, double seconds)
{
    cout << section << "," << test << "," << dip << "," << aspect
	 << "," << component << "," << value << "," << reference
	 << "," << error << "," << nphi << "," << nq << "," << ntheta
	 << "," << nk << "," << na << "," << seconds << endl;
}

int main(int argc, char** argv)
{
    double c11=47.0, c12=8.0, c13=5.0, c33=34.0, c44=17.0;
    double factor_w=0.4;
    int nphi=360, nq=64, ntheta=400;
    bool fourier_only=false, gotError=false;
    string progname("benchmark_kanaun09");

    for (int i=1; i<argc; i++){
	if (argv[i][0]=='-'){
	    switch(argv[i][1]){
	    case 'C':
		if (sscanf(&argv[i][2],"%lg/%lg/%lg/%lg/%lg",
			   &c11,&c12,&c13,&c33,&c44)!=5){
		    cerr << progname << ": incomplete -C option\n";
		    gotError=true;
		}
		break;
	    case 'P':
		nphi=atoi(&argv[i][2]);
		break;
	    case 'Q':
		nq=atoi(&argv[i][2]);
		break;
	    case 'N':
		ntheta=atoi(&argv[i][2]);
		break;
	    case 'W':
		factor_w=atof(&argv[i][2]);
		break;
	    case 'F':
		fourier_only=true;
		break;
	    default:
		cerr << progname << ": unknown option detected: " << argv[i] << '\n';
		gotError=true;
		break;
	    }
	}else{
	    cerr << progname << ": unknown option detected: " << argv[i] << '\n';
	    gotError=true;
	}
    }
    if (nphi<4 || nq<4 || ntheta<4) gotError=true;
    if (gotError){
	cerr << progname << ": incorrect usage detected.\n";
	return 1;
    }

    Array1d<double> c(5);
    c(1)=c11; c(2)=c12; c(3)=c13; c(4)=c33; c(5)=c44;
    Tensor4 C(c);

    cout.precision(12);
    cout << "section,test,dip_deg,aspect,component,value,reference,"
	 << "relative_error,nphi,nq,ntheta,kanaun_points,adaptive_points,seconds\n";

    // Eq. (82): circular crack in an isotropic medium.
    Array1d<double> ci(2);
    ci(1)=80.0; ci(2)=30.0;
    Tensor4 Ciso(ci);
    auto start=chrono::steady_clock::now();
    Kanaun09 Kiso(Ciso,1.0,1.0,0.0,0.0,0.0,nphi,nq);
    auto stop=chrono::steady_clock::now();
    double seconds=chrono::duration<double>(stop-start).count();
    Array2d<double> T0(3,3), Tref(3,3);
    Kiso.setT0(T0);
    kanaun09_isotropic_T0(Ciso,1.0,Tref);
    for (int i=1; i<=3; i++){
	double error=abs(T0(i,i)-Tref(i,i))/abs(Tref(i,i));
	string comp=(i==1) ? "T11" : ((i==2) ? "T22" : "T33");
	print_row("analytic","isotropic",0.0,0.0,comp,T0(i,i),Tref(i,i),
		  error,nphi,nq,ntheta,Kiso.numIntegPoints(),0,seconds);
    }

    // Eqs. (76)-(77): circular crack in the isotropy plane of a TI host.
    start=chrono::steady_clock::now();
    Kanaun09 Kti(C,1.0,1.0,0.0,0.0,0.0,nphi,nq);
    stop=chrono::steady_clock::now();
    seconds=chrono::duration<double>(stop-start).count();
    Kti.setT0(T0);
    kanaun09_transverse_T0(C,1.0,Tref);
    for (int i=1; i<=3; i++){
	double error=abs(T0(i,i)-Tref(i,i))/abs(Tref(i,i));
	string comp=(i==1) ? "T11" : ((i==2) ? "T22" : "T33");
	print_row("analytic","transverse",0.0,0.0,comp,T0(i,i),Tref(i,i),
		  error,nphi,nq,ntheta,Kti.numIntegPoints(),0,seconds);
    }

    // Numerical example of Section 4.1.3 and Fig. 4 of Kanaun and Levin.
    Array1d<double> cp(5);
    cp(1)=50.0; cp(2)=10.0; cp(3)=30.0; cp(4)=100.0; cp(5)=60.0;
    Tensor4 Cpaper(cp);
    for (int idip=0; idip<=36; idip++){
	double dip=5.0*idip;
	start=chrono::steady_clock::now();
	// Kanaun's positive theta in Eq. (20) has the opposite sign from
	// the parent's positive intrinsic roll angle gamma.
	Kanaun09 Kpaper(Cpaper,1.0,1.0,0.0,0.0,-dip,nphi,nq);
	stop=chrono::steady_clock::now();
	seconds=chrono::duration<double>(stop-start).count();
	Array2d<double> Mp(3,3);
	Kpaper.setM(Mp);
	print_row("paper_example","figure4",dip,0.0,"M11",Mp(1,1),0.0,0.0,
		  nphi,nq,ntheta,Kpaper.numIntegPoints(),0,seconds);
	print_row("paper_example","figure4",dip,0.0,"M22",Mp(2,2),0.0,0.0,
		  nphi,nq,ntheta,Kpaper.numIntegPoints(),0,seconds);
	print_row("paper_example","figure4",dip,0.0,"M33",Mp(3,3),0.0,0.0,
		  nphi,nq,ntheta,Kpaper.numIntegPoints(),0,seconds);
	print_row("paper_example","figure4",dip,0.0,"M23",Mp(2,3),0.0,0.0,
		  nphi,nq,ntheta,Kpaper.numIntegPoints(),0,seconds);
    }

    if (!fourier_only){
	// The 2009 solution has zero thickness. The comparable finite-spheroid
	// quantity is aspect_ratio*H, not H itself.
	double aspects[5]={1e-1,3e-2,1e-2,3e-3,1e-3};
	double dips[4]={0.0,30.0,60.0,90.0};
	for (int idip=0; idip<4; idip++){
	    double dip=dips[idip];
	    start=chrono::steady_clock::now();
	    Kanaun09 K(C,1.0,1.0,0.0,0.0,dip,nphi,nq);
	    stop=chrono::steady_clock::now();
	    double kanaun_seconds=chrono::duration<double>(stop-start).count();
	    Tensor4 H0;
	    K.setCrackCompliance(H0);

	    for (int ia=0; ia<5; ia++){
		double aspect=aspects[ia];
		Ellipsoid E(1.0,1.0,aspect,0.0,0.0,dip);
		start=chrono::steady_clock::now();
		HillTensor Hill(C,E,ntheta,2*ntheta,factor_w,false,20);
		Tensor4 P;
		Hill.setTensor(P);
		Tensor4 H=emt::dry_compliance_contribution(C,P);
		H *= aspect;
		stop=chrono::steady_clock::now();
		seconds=kanaun_seconds+chrono::duration<double>(stop-start).count();
		double rms=relative_rms(H,H0);
		double maxdiff=max_abs_diff(H,H0);
		print_row("finite_limit","adaptive",dip,aspect,"relative_rms",
			  rms,0.0,0.0,nphi,nq,ntheta,K.numIntegPoints(),
			  Hill.numIntegPoints(),seconds);
		print_row("finite_limit","adaptive",dip,aspect,"max_abs",
			  maxdiff,0.0,0.0,nphi,nq,ntheta,K.numIntegPoints(),
			  Hill.numIntegPoints(),seconds);
	    }
	}
    }

    return 0;
}
