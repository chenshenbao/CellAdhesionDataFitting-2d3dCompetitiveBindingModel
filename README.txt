1: dpn.m --- defines the following differential equations of 2D molecular binding and dissociation kinetics: 

   dp(n)/dt = ac*mr*ml*kf*p(n-1) - (ac*mr*ml*kf + n*kr)*p(n) + (n+1)*kr*p(n+1)

   dpnm.m ---  defines the following differential equations of 2D and 3D competitive binding and dissociation kinetics: 

   dp(n,m)/dt = kf2*ml*[ac*mr-m-(n-1)]*p(n-1,m) + kr2*(n+1).*p(n+1,m) + kf3*c*[ac*mr-(m-1)-n]*p(n,m-1) + kr3*(m+1).*p(n,m+1) - [kf2*ml*(ac*mr-m-n)+kr2*n + kf3*c*(ac*mr-m-n)+kr3*m]*p(n,m);

   equations on and outside the boundaries need to be modified. We find the boundaries and adjust the equations by coefficient matrices.

2: Solving of the equations is based on the MATLAB ode45 solver;

3: Based on the p matrix given by ode45, we can get cell adhesion probabilities by excluding p(m,0);

4: Fitting progress is based on the MATLAB  function lsqnonlin();

5: run main.m to fit experimental data.
