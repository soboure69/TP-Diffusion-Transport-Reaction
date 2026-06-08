function f = fef1l_x(xl,xr)
%----------------------------------------------------------
% element vector for the source term f(x)=x on [xl,xr]
%
% f_i = int_{xl}^{xr} x * N_i(x) dx
% for linear shape functions N1, N2
%----------------------------------------------------------

eleng = xr - xl;

f = (eleng/6) * [2*xl + xr;
                 xl + 2*xr];
end