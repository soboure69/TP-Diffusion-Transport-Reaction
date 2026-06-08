function [eL2, eH1semi, emax] = compute_error(xnodes, uh, uexact, duexact)

%-------------------------------------------------------------------
%  Compute error norms between the P1 FE solution and the exact one.
%
%  Inputs:
%    xnodes  - nodal coordinates (n+1 x 1)
%    uh      - FE nodal values (n+1 x 1)
%    uexact  - function handle u_exact(x)
%    duexact - (optional) function handle u_exact'(x) for the H1 seminorm
%
%  Outputs:
%    eL2     - L2 norm of the error    sqrt( int (u-uh)^2 )
%    eH1semi - H1 seminorm of error    sqrt( int (u'-uh')^2 )  (NaN if duexact missing)
%    emax    - max nodal error         max |u(x_i) - uh_i|
%-------------------------------------------------------------------

if nargin < 4, duexact = []; end

n = length(xnodes) - 1;

% 3-point Gauss on [-1,1] (exact up to degree 5)
gp = [-sqrt(3/5), 0, sqrt(3/5)];
gw = [5/9, 8/9, 5/9];

eL2sq = 0;
eH1sq = 0;
for iel = 1:n
  xl = xnodes(iel);
  xr = xnodes(iel+1);
  L  = xr - xl;
  uhl = uh(iel);
  uhr = uh(iel+1);
  duh = (uhr - uhl)/L;            % FE derivative (constant per element)
  for g = 1:3
    xi  = gp(g);
    N1  = 0.5*(1-xi);
    N2  = 0.5*(1+xi);
    x   = N1*xl + N2*xr;
    jac = L/2;
    uh_g = N1*uhl + N2*uhr;
    eL2sq = eL2sq + gw(g) * (uexact(x) - uh_g)^2 * jac;
    if ~isempty(duexact)
      eH1sq = eH1sq + gw(g) * (duexact(x) - duh)^2 * jac;
    end
  end
end

eL2 = sqrt(eL2sq);
if ~isempty(duexact)
  eH1semi = sqrt(eH1sq);
else
  eH1semi = NaN;
end

% max nodal error
uex_nodes = arrayfun(uexact, xnodes);
emax = max(abs(uex_nodes - uh));
