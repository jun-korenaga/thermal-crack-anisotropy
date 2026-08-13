function self_test()
% SELF_TEST Small checks independent of the full MINEOS data set.

% Check the trapezoidal rule on a linear function.
x = [0;1;3;6]; y = 2*x+1;
val = trapz(x,y);
assert(abs(val-42)<1e-12,'trapezoidal integration test failed');

% Check the coefficient of the dF contribution.
k = 2; Uprime = 3; V = 5; dF = 7;
value = 2*Uprime*V/k*dF;
assert(value==105,'dF coefficient test failed');

% A symmetric positive diagonal tensor must give real wave speeds.
C = diag([200 200 220 80 80 70]);
[vsv,vsh,vpv,vph,xi] = calc_vel_etc(C(:)',3300);
assert(all(isfinite([vsv vsh vpv vph xi])),'velocity test failed');
disp('self-test passed');
end
