%% Q5a - Transport dominant (sigma=0, b>>mu) : Galerkin standard
%  -mu u'' + b u' = 0, u(0)=0, u(1)=1   avec  mu=0.01, b=1
%  Mise en application DIRECTE des routines fournies (sans solve_dtr).
%
%  Forme du template fourni :  a u'' + b u' + c u = f
%  Ici : a=acoef=-mu=-0.01, b=bcoef=1, c=ccoef=0, f=0.
%  But : montrer les oscillations de Galerkin quand Pe_h = b*h/(2*mu) > 1.

clear; close all;
if exist('OCTAVE_VERSION','builtin'), graphics_toolkit('fltk'); warning('off','all'); end

%--- donnees physiques ---
mu = 0.01; b_coef = 1; sigma = 0;
acoef = -mu;  bcoef = b_coef;  ccoef = sigma;   % coefficients du template

nnel = 2; ndof = 1;
n_values = [10, 40, 60];

%--- solution exacte (analytique) en reference ---
x_exact = linspace(0,1,500)';
alpha = b_coef/mu;
u_exact = (exp(alpha*x_exact) - 1) / (exp(alpha) - 1);

figure();
plot(x_exact, u_exact, 'k-', 'LineWidth', 2); hold on;
colors = {'r', 'b', 'g'};
markers = {'o', 's', 'd'};

for i = 1:length(n_values)
  nel = n_values(i);
  nnode = nel + 1;  sdof = nnode*ndof;
  hh = 1/nel;

  %--- coordonnees et connectivite ---
  gcoord = zeros(nnode,1);
  for k1 = 1:nnode, gcoord(k1) = (k1-1)*hh; end
  nodes = zeros(nel,2);
  for iel = 1:nel, nodes(iel,1)=iel; nodes(iel,2)=iel+1; end

  %--- conditions aux limites ---
  bcdof = [1, nnode];  bcval = [0, 1];   % u(0)=0, u(1)=1

  %--- assemblage avec les routines fournies (f=0) ---
  ff = zeros(sdof,1);
  kk = zeros(sdof,sdof);
  for iel = 1:nel
    xl = gcoord(nodes(iel,1));  xr = gcoord(nodes(iel,2));
    eleng = xr - xl;
    index = feeldof1(iel,nnel,ndof);          % fourni
    ke = feode2l(acoef,bcoef,ccoef,eleng);    % fourni
    fe = [0; 0];                              % f = 0
    [kk,ff] = feasmbl2(kk,ff,ke,fe,index);    % fourni
  end
  [kk,ff] = feaplyc2(kk,ff,bcdof,bcval);      % fourni
  uh = kk \ ff;

  plot(gcoord, uh, [colors{i} markers{i} '-'], 'MarkerSize', 4, 'LineWidth', 1);
  Pe = b_coef*hh/(2*mu);
  fprintf('n=%d, h=%.4f, Pe_h=%.2f, max(uh)=%.4f, min(uh)=%.4f\n', ...
          nel, hh, Pe, max(uh), min(uh));
end

xlabel('x'); ylabel('u(x)');
title(sprintf('Q5a: Transport dominant - Galerkin (\\mu=%.2f, b=%g)', mu, b_coef));
legend('Exacte', sprintf('n=%d',n_values(1)), sprintf('n=%d',n_values(2)), ...
       sprintf('n=%d',n_values(3)), 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q5a_transport_galerkin.png');
fprintf('\nQ5a done.\n');
fprintf('Oscillations appear when the element Peclet number Pe_h = b*h/(2*mu) > 1.\n');
fprintf('For b=%g, mu=%g: h_crit = 2*mu/b = %g (n_crit = %g).\n', ...
        b_coef, mu, 2*mu/b_coef, ceil(b_coef/(2*mu)));
