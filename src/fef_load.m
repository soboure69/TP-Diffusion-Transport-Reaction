function [f]=fef_load(fun,xl,xr)

%-------------------------------------------------------------------
%  Purpose:
%     Element load vector  f_i = \int_{xl}^{xr} fun(x) N_i(x) dx
%     for a linear (P1) element, using 2-point Gauss quadrature
%     (exact for polynomial integrands up to degree 3).
%
%  Synopsis:
%     [f]=fef_load(fun,xl,xr)
%
%  Variable Description:
%     f   - element vector (size 2x1)
%     fun - function handle of the source term f(x)
%     xl  - coordinate of the left node
%     xr  - coordinate of the right node
%-------------------------------------------------------------------

 eleng = xr-xl;                 % element length
 % 2-point Gauss points / weights on the reference interval [-1,1]
 gp = [-1/sqrt(3), 1/sqrt(3)];
 gw = [1, 1];

 f = [0; 0];
 for g=1:2
   xi = gp(g);
   % shape functions in terms of the reference coordinate xi in [-1,1]
   N1 = 0.5*(1-xi);
   N2 = 0.5*(1+xi);
   x  = N1*xl + N2*xr;          % physical coordinate
   jac = eleng/2;               % Jacobian dx/dxi
   fx = fun(x);
   f(1) = f(1) + gw(g)*fx*N1*jac;
   f(2) = f(2) + gw(g)*fx*N2*jac;
 end
