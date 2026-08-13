/*
 * kanaun09.cc
 *
 * Summer 2026
 * Jun Korenaga
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "kanaun09.h"

using namespace std;

Kanaun09::Kanaun09(const Tensor4& _C, double _a1, double _a2,
		   double _alpha, double _beta, double _gamma,
		   int _nphi, int _nq)
    : C(_C), a1(_a1), a2(_a2),
      alpha(_alpha), beta(_beta), gamma(_gamma),
      nphi(_nphi), nq(_nq), n_int_points(0)
{
    if (a1<=0.0 || a2<=0.0 || nphi<4 || nq<4){
	cerr << "Kanaun09: invalid input detected.\n";
	exit(1);
    }

    R.resize(3,3);
    T0.resize(3,3);
    M.resize(3,3);
    xq.resize(nq);
    wq.resize(nq);

    set_R();
    set_C_local();
    set_quadrature();
    set_T0();
    set_M();
    set_H();
}

void Kanaun09::set_R()
{
    double deg2rad=M_PI/180.0;
    alpha *= deg2rad;
    beta *= deg2rad;
    gamma *= deg2rad;

    double ca=cos(alpha), sa=sin(alpha);
    double cb=cos(beta), sb=sin(beta);
    double cg=cos(gamma), sg=sin(gamma);

    R(1,1) = ca*cb;
    R(1,2) = ca*sb*sg-sa*cg;
    R(1,3) = ca*sb*cg+sa*sg;
    R(2,1) = sa*cb;
    R(2,2) = sa*sb*sg+ca*cg;
    R(2,3) = sa*sb*cg-ca*sg;
    R(3,1) = -sb;
    R(3,2) = cb*sg;
    R(3,3) = cb*cg;
}

void Kanaun09::set_C_local()
{
    for (int i=1; i<=3; i++){
	for (int j=1; j<=3; j++){
	    for (int k=1; k<=3; k++){
		for (int l=1; l<=3; l++){
		    double val=0.0;
		    for (int a=1; a<=3; a++){
			for (int b=1; b<=3; b++){
			    for (int c=1; c<=3; c++){
				for (int d=1; d<=3; d++){
				    val += R(a,i)*R(b,j)*R(c,k)*R(d,l)
					*C.val(a,b,c,d);
				}
			    }
			}
		    }
		    C_local.set(i,j,k,l,val);
		}
	    }
	}
    }
}

void Kanaun09::set_quadrature()
{
    // Gauss-Legendre abscissae and weights on [-1,1].
    int m=(nq+1)/2;
    double eps=1e-15;

    for (int i=1; i<=m; i++){
	double z=cos(M_PI*(i-0.25)/(nq+0.5));
	double z1, p1, p2, pp;
	do {
	    p1=1.0;
	    p2=0.0;
	    for (int j=1; j<=nq; j++){
		double p3=p2;
		p2=p1;
		p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	    }
	    pp=nq*(z*p1-p2)/(z*z-1.0);
	    z1=z;
	    z=z1-p1/pp;
	} while (abs(z-z1)>eps);

	xq(i)=-z;
	xq(nq+1-i)=z;
	wq(i)=2.0/((1.0-z*z)*pp*pp);
	wq(nq+1-i)=wq(i);
    }
}

void Kanaun09::invert3(const Array2d<double>& A,
		       Array2d<double>& invA) const
{
    Array2d<double> tmp(3,3);
    tmp=A;
    invA.resize(3,3);
    invA=0.0;
    for (int i=1; i<=3; i++) invA(i,i)=1.0;

    for (int i=1; i<=3; i++){
	int pivot=i;
	for (int j=i+1; j<=3; j++){
	    if (abs(tmp(j,i))>abs(tmp(pivot,i))) pivot=j;
	}
	if (abs(tmp(pivot,i))<1e-30){
	    cerr << "Kanaun09::invert3 - singular matrix detected.\n";
	    exit(1);
	}
	if (pivot!=i){
	    for (int j=1; j<=3; j++){
		swap(tmp(i,j),tmp(pivot,j));
		swap(invA(i,j),invA(pivot,j));
	    }
	}

	double scale=1.0/tmp(i,i);
	for (int j=1; j<=3; j++){
	    tmp(i,j) *= scale;
	    invA(i,j) *= scale;
	}
	for (int j=1; j<=3; j++){
	    if (j!=i){
		double factor=tmp(j,i);
		for (int k=1; k<=3; k++){
		    tmp(j,k) -= factor*tmp(i,k);
		    invA(j,k) -= factor*invA(i,k);
		}
	    }
	}
    }
}

void Kanaun09::calc_S(double q1, double q2, double q3,
		      Array2d<double>& S) const
{
    Array1d<double> q(3);
    q(1)=q1; q(2)=q2; q(3)=q3;

    Array2d<double> L(3,3), G(3,3);
    for (int i=1; i<=3; i++){
	for (int j=1; j<=3; j++){
	    double val=0.0;
	    for (int l=1; l<=3; l++){
		for (int m=1; m<=3; m++){
		    val += q(l)*C_local.val(l,i,j,m)*q(m);
		}
	    }
	    L(i,j)=val;
	}
    }
    invert3(L,G);

    S.resize(3,3);
    for (int i=1; i<=3; i++){
	for (int j=1; j<=3; j++){
	    double val=0.0;
	    for (int n=1; n<=3; n++){
		double left=0.0;
		for (int m=1; m<=3; m++){
		    left += C_local.val(3,i,m,n)*q(m);
		}
		for (int p=1; p<=3; p++){
		    double right=0.0;
		    for (int r=1; r<=3; r++){
			right += C_local.val(p,r,j,3)*q(r);
		    }
		    val += left*G(n,p)*right;
		}
	    }
	    S(i,j)=val-C_local.val(3,i,j,3);
	}
    }
}

void Kanaun09::set_T0()
{
    // Eqs. (23) and (32). The tangent map transforms q3=(-inf,inf)
    // to the Gauss-Legendre interval x=(-1,1).
    T0=0.0;
    double dphi=2.0*M_PI/nphi;
    Array2d<double> S(3,3);

    for (int iphi=1; iphi<=nphi; iphi++){
	double phi=(iphi-0.5)*dphi;
	double q1=cos(phi)/a1;
	double q2=sin(phi)/a2;
	double qbar=sqrt(q1*q1+q2*q2);
	Array2d<double> sum(3,3);
	sum=0.0;

	for (int iq=1; iq<=nq; iq++){
	    double angle=0.5*M_PI*xq(iq);
	    double cosine=cos(angle);
	    double q3=qbar*tan(angle);
	    double jac=0.5*M_PI*qbar/(cosine*cosine);
	    calc_S(q1,q2,q3,S);
	    for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++){
		    sum(i,j) += wq(iq)*jac*S(i,j);
		}
	    }
	    n_int_points++;
	}
	for (int i=1; i<=3; i++){
	    for (int j=1; j<=3; j++){
		T0(i,j) -= dphi*sum(i,j)/(8.0*M_PI);
	    }
	}
    }
}

void Kanaun09::set_M()
{
    invert3(T0,M);
    M /= sqrt(M_PI*a1*a2);
}

void Kanaun09::set_H()
{
    // Eq. (132), converted from crack number density to the limit of
    // aspect_ratio * H, where porosity=(4*pi/3)*n*a^3*aspect_ratio.
    Tensor4 H_local;
    H_local=0.0;
    double fac=sqrt(M_PI)/2.0;

    for (int i=1; i<=3; i++){
	for (int j=1; j<=3; j++){
	    for (int k=1; k<=3; k++){
		for (int l=1; l<=3; l++){
		    double mi=(i==3) ? 1.0 : 0.0;
		    double mj=(j==3) ? 1.0 : 0.0;
		    double mk=(k==3) ? 1.0 : 0.0;
		    double ml=(l==3) ? 1.0 : 0.0;
		    double val=0.25*(mi*M(j,k)*ml+mj*M(i,k)*ml
				     +mi*M(j,l)*mk+mj*M(i,l)*mk);
		    H_local.set(i,j,k,l,fac*val);
		}
	    }
	}
    }

    H=0.0;
    for (int i=1; i<=3; i++){
	for (int j=1; j<=3; j++){
	    for (int k=1; k<=3; k++){
		for (int l=1; l<=3; l++){
		    double val=0.0;
		    for (int a=1; a<=3; a++){
			for (int b=1; b<=3; b++){
			    for (int c=1; c<=3; c++){
				for (int d=1; d<=3; d++){
				    val += R(i,a)*R(j,b)*R(k,c)*R(l,d)
					*H_local.val(a,b,c,d);
				}
			    }
			}
		    }
		    H.set(i,j,k,l,val);
		}
	    }
	}
    }
}

void kanaun09_isotropic_T0(const Tensor4& C, double a,
			   Array2d<double>& T0)
{
    double c11=C.val(1,1), c12=C.val(1,2), mu=C.val(4,4);
    double nu=c12/(c11+c12);
    T0.resize(3,3);
    T0=0.0;
    T0(1,1)=T0(2,2)=M_PI*mu*(2.0-nu)/(8.0*a*(1.0-nu));
    T0(3,3)=M_PI*mu/(4.0*a*(1.0-nu));
}

void kanaun09_transverse_T0(const Tensor4& C, double a,
			    Array2d<double>& T0)
{
    // Eqs. (57), (61), (76), and (77), for a circular crack in the
    // isotropy plane.
    double kappa=0.5*(C.val(1,1)+C.val(1,2));
    double m=0.5*(C.val(1,1)-C.val(1,2));
    double n=C.val(3,3), ell=C.val(1,3), mu=C.val(4,4);
    double aa=(n*(kappa+m)-2.0*ell*mu-ell*ell)/(mu*n);
    double bb=(kappa+m)/n;
    double disc=aa*aa-4.0*bb;
    if (disc<0.0){
	cerr << "kanaun09_transverse_T0: complex roots are not supported.\n";
	exit(1);
    }
    double xi2=0.5*(aa+sqrt(disc));
    double xi3=0.5*(aa-sqrt(disc));
    double den=n*(sqrt(xi2)+sqrt(xi3));
    double term=(n*(kappa+m)-ell*ell)/den;

    T0.resize(3,3);
    T0=0.0;
    T0(1,1)=T0(2,2)=M_PI*(sqrt(mu*m)+term)/(8.0*a);
    T0(3,3)=M_PI*(n*(kappa+m)-ell*ell)
	/(4.0*a*n*sqrt(xi2*xi3)*(sqrt(xi2)+sqrt(xi3)));
}
