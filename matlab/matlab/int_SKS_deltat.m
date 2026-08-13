function val = int_SKS_deltat(i,rho,C)
% INT_SKS_DELTAT Integrand for maximum SKS delay time.

Ctmp = reshape(C(i,:,:),[6,6]);
L = 0.5*(Ctmp(4,4)+Ctmp(5,5));
Gc = 0.5*(Ctmp(5,5)-Ctmp(4,4));

if L>0
  val = sqrt(rho(i)/L)*(Gc/L);
else
  val = 0;
end
end
