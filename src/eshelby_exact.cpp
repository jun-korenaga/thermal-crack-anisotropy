/*
 * eshelby_exact.cc
 *
 * Jun Korenaga
 * Summer 2024
 */

#include "emt/hill_tensor.hpp"

using namespace std;

Tensor4 eshelby_isotropic_exact(const Tensor4& C, double a1, double a3)
{
    Tensor4 s;
    s=0.0;

    double c11, c12, nu, gamma, g, nu_1;

    c11 = C.val(1,1); c12 = C.val(1,2);
    nu = c12/(c11+c12);
    gamma = a3/a1;
    nu_1 = 1.0-nu;

    double val;

    if (gamma==1.0){
	val = (7.0-5.0*nu)/(15.0*nu_1);
	s.set(1,1,1,1,val);
	s.set(2,2,2,2,val);
	s.set(3,3,3,3,val);

	val = (4.0-5.0*nu)/(15.0*nu_1);
	s.set(1,2,1,2,val);
	s.set(2,3,2,3,val);
	s.set(3,1,3,1,val);

	val = (5.0*nu-1.0)/(15.0*nu_1);
	s.set(1,1,2,2,val);
	s.set(2,2,3,3,val);
	s.set(3,3,1,1,val);
	s.set(1,1,3,3,val);
	s.set(2,2,1,1,val);
	s.set(3,3,2,2,val);
    }else{
	double gamma2 = gamma*gamma;
	if (gamma<1.0){
	    g = gamma/pow(1.0-gamma2,1.5)
		*(acos(gamma)-gamma*sqrt(1.0-gamma2));
       }else{
	    g = gamma/pow(gamma2-1.0,1.5)
		*(gamma*sqrt(gamma2-1.0)-acosh(gamma));
	}
	double gamma2_1 = 1.0-gamma2;

	val = -3.0*gamma2/(8.0*nu_1*gamma2_1)
	    + 1.0/(4.0*nu_1)*(1.0-2.0*nu+9.0/(4.0*gamma2_1))*g;
	s.set(1,1,1,1,val);
	s.set(2,2,2,2,val);

	val = 1.0/nu_1*(2.0-nu-1.0/gamma2_1)
	    + 1.0/(2.0*nu_1)*(-4.0+2.0*nu+3.0/gamma2_1)*g;
	s.set(3,3,3,3,val);

	val = 1.0/(8.0*nu_1)*(1.0-1.0/gamma2_1)
	    + 1.0/(16.0*nu_1)*(-4.0*(1.0-2.0*nu)+3.0/gamma2_1)*g;
	s.set(1,1,2,2,val);
	s.set(2,2,1,1,val);

	val = gamma2/(2.0*nu_1*gamma2_1)
	    -1.0/(4.0*nu_1)*(1.0-2.0*nu+3.0*gamma2/gamma2_1)*g;
	s.set(1,1,3,3,val);
	s.set(2,2,3,3,val);

	val = 1.0/(2.0*nu_1)*(-(1.0-2.0*nu)+1.0/gamma2_1)
	    + 1.0/(4.0*nu_1)*(2.0*(1.0-2.0*nu)-3.0/gamma2_1)*g;
	s.set(3,3,1,1,val);
	s.set(3,3,2,2,val);

	val = -gamma2/(8.0*nu_1*gamma2_1)
	    + 1.0/(16.0*nu_1)*(4.0*(1.0-2.0*nu)+3.0/gamma2_1)*g;
	s.set(1,2,1,2,val);

	val = 1.0/(4.0*nu_1)*(1.0-2.0*nu+(1.0+gamma2)/gamma2_1)
	    - 1.0/(8.0*nu_1)*(1.0-2.0*nu+3.0*(1.0+gamma2)/gamma2_1)*g;
	s.set(1,3,1,3,val);
	s.set(2,3,2,3,val);
    }
    
    return s;
}

Tensor4 eshelby_transverse_exact(const Tensor4& C, double a1, double a3)
{
    Tensor4 s;
    s=0.0;

    double c11, c13, c33, c44, c66, a, c, C1, D, G, F;
    Array1d<double> k(2), A(2), B(5), I1(3), I2(3);

    c11 = C.val(1,1); c13 = C.val(1,3); 
    c33 = C.val(3,3); c44 = C.val(4,4); c66 = C.val(6,6);
    a = a1; c = a3; 

    B(3) = sqrt(c66/c44);
    C1 = sqrt(c11*c33);
    B(4) = sqrt((C1-c13)*(C1+c13+2.0*c44)/(c33*c44));
    B(5) = sqrt((C1+c13)*(C1-c13-2.0*c44)/(c33*c44));
    B(1) = 0.5*(B(4)+B(5));
    B(2) = 0.5*(B(4)-B(5));
    D = 1.0/(4.0*M_PI*c44*B(3));

    for (int i=1; i<=2; i++){
	double Bi2 = B(i)*B(i);
	k(i) = (c11/Bi2-c44)/(c13+c44);
	double fac = (i==1) ? -1.0 : 1.0;
	A(i) = fac*(c44-c33*Bi2)/(8.0*M_PI*c33*c44*(B(1)*B(1)-B(2)*B(2))*Bi2);
    }
    for (int i=1; i<=3; i++){
	if (B(i)*c>a){
	    G = sqrt(B(i)*B(i)*c*c-a*a);
	    F = acosh(B(i)*c/a);
	    double G3 = pow(G,3.0);
	    I1(i) = (2.0*M_PI*c/G3)*(B(i)*c*G-a*a*F);
	    I2(i) = (4.0*M_PI*a*a*c/G3)*(F-G/(B(i)*c));
	}else{
	    G = sqrt(a*a-B(i)*B(i)*c*c);
	    F = acos(B(i)*c/a);
	    double G3 = pow(G,3.0);
	    I1(i) = -(2.0*M_PI*c/G3)*(B(i)*c*G-a*a*F);
	    I2(i) = -(4.0*M_PI*a*a*c/G3)*(F-G/(B(i)*c));
	}
    }

    double sum1, sum2, val;

    // s_1111, s_1122
    sum1=0.0; sum2=0.0;
    for (int i=1; i<=2; i++){
	sum1 += c44*(1+k(i))*A(i)*pow(B(i),3.0)*I1(i);
	sum2 += A(i)*B(i)*I1(i);
    }
    val = 2.0*sum1 - c66*sum2 + 0.5*D*c66*I1(3);
    s.set(1,1,1,1,val); // (1,1)
    s.set(2,2,2,2,val); // (2,2)
    val = 2.0*sum1 - 3.0*c66*sum2 - 0.5*D*c66*I1(3);
    s.set(1,1,2,2,val); // (1,2)
    s.set(2,2,1,1,val); // (2,1)

    // s_3333
    sum1=0.0;
    for (int i=1; i<=2; i++){
	double Bi2 = B(i)*B(i);
	double Bi3 = Bi2*B(i);
	sum1 += Bi3*k(i)*A(i)*(c13-c33*k(i)*Bi2)*I2(i);
    }
    val = 2.0*sum1;
    s.set(3,3,3,3,val); // (3,3)

    // s_1133
    sum1=0.0;
    for (int i=1; i<=2; i++){
	sum1 += B(i)*A(i)*(c13-c33*k(i)*B(i)*B(i))*I1(i);
    }
    val = 2.0*sum1;
    s.set(1,1,3,3,val); // (1,3)
    s.set(2,2,3,3,val); // (2,3)

    // s_3311
    sum1=0.0; sum2=0.0;
    for (int i=1; i<=2; i++){
	double Bi3 = pow(B(i),3.0);
	double Bi5 = pow(B(i),5.0);
	sum1 += c44*Bi5*k(i)*A(i)*(1.0+k(i))*I2(i);
	sum2 += c66*Bi3*k(i)*A(i)*I2(i);
    }
    val = 2.0*sum1 - 2.0*sum2;
    s.set(3,3,1,1,val); // (3,1)
    s.set(3,3,2,2,val); // (3,2)

    // s_1212
    sum1=0.0;
    for (int i=1; i<=2; i++){
	sum1 += A(i)*B(i)*I1(i);
    }
    val = c66*sum1 + 0.5*D*c66*I1(3);
    s.set(1,2,1,2,val); // (6,6)

    // s_1313
    sum1=0.0;
    for (int i=1; i<=2; i++){
	sum1 += A(i)*pow(B(i),3.0)*(1.0+k(i))*(I2(i)-2.0*k(i)*I1(i));
    }
    val = 0.5*c44*sum1 + 0.25*D*c44*I2(3)*B(3)*B(3);
    s.set(1,3,1,3,val); // (5,5)
    s.set(2,3,2,3,val); // (4,4)

    return s;
}
