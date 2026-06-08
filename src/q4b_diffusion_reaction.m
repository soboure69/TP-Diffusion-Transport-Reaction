%% Q4b - Diffusion-Reaction : -u'' + u = x, u(0)=0, u(1)=1
%  Mise en application DIRECTE des routines fournies (comme Diff_Trans_Reac1D.m)
%  pour les parametres specifiques de Q4b, puis trace de la solution approchee u_h.
%
%  Forme du template fourni :  a u'' + b u' + c u = f
%  Ici : a=acoef=-1, b=bcoef=0, c=ccoef=1, f(x)=x  =>  -u'' + u = x
%
%  Remarque : f(x)=x n'est PAS constante, donc fef1l (qui ne fait que f=1)
%  ne convient pas. On assemble le vecteur charge elementaire a la main par
%  quadrature de Gauss a 2 points (exacte ici). L'assemblage de la matrice
%  et les conditions aux limites utilisent les routines fournies.

clear; close all;
if exist('OCTAVE_VERSION','builtin'), graphics_toolkit('fltk'); warning('off','all'); end

%--- parametres de controle ---
nel   = 20;            % nombre d'elements
nnel  = 2;             % noeuds par element
ndof  = 1;             % ddl par noeud
nnode = nel + 1;       % nombre total de noeuds
sdof  = nnode*ndof;    % ddl total

%--- coordonnees nodales ---
h = 1/nel;
for i = 1:nnode
  gcoord(i) = (i-1)*h;
end

%--- connectivite ---
for iel = 1:nel
  nodes(iel,1) = iel;
  nodes(iel,2) = iel+1;
end

%--- coefficients de l'EDO ---
acoef = -1;   % a (= -mu)
bcoef =  0;   % b
ccoef =  1;   % c (= sigma)
ffun  = @(x) x;   % second membre f(x)=x

%--- conditions aux limites ---
bcdof(1) = 1;     bcval(1) = 0;     % u(0) = 0
bcdof(2) = nnode; bcval(2) = 1;     % u(1) = 1

%--- initialisation ---
ff = zeros(sdof,1);
kk = zeros(sdof,sdof);

%--- points de Gauss a 2 points sur [-1,1] ---
gp = [-1/sqrt(3), 1/sqrt(3)];
gw = [1, 1];

%--- assemblage ---
for iel = 1:nel
  nl = nodes(iel,1); nr = nodes(iel,2);
  xl = gcoord(nl);   xr = gcoord(nr);
  eleng = xr - xl;
  index = feeldof1(iel,nnel,ndof);          % ddl de l'element (fourni)
  k = feode2l(acoef,bcoef,ccoef,eleng);     % matrice elementaire (fourni)

  % vecteur charge elementaire pour f(x) quelconque (Gauss 2 pts, exact ici)
  f = [0; 0];
  for g = 1:2
    xi  = gp(g);
    N1  = 0.5*(1-xi);
    N2  = 0.5*(1+xi);
    xg  = N1*xl + N2*xr;
    jac = eleng/2;
    fx  = ffun(xg);
    f(1) = f(1) + gw(g)*fx*N1*jac;
    f(2) = f(2) + gw(g)*fx*N2*jac;
  end

  [kk,ff] = feasmbl2(kk,ff,k,f,index);      % assemblage (fourni)
end

%--- application des conditions de Dirichlet (fourni) ---
[kk,ff] = feaplyc2(kk,ff,bcdof,bcval);

%--- resolution ---
fsol = kk \ ff;

%--- representation graphique de la solution approchee u_h ---
figure();
plot(gcoord, fsol, 'ro-', 'MarkerSize', 5, 'LineWidth', 1.2);
xlabel('x'); ylabel('u_h(x)');
title(sprintf('Q4b: solution approchee u_h (-u'''' + u = x, n=%d)', nel));
legend('Solution EF u_h (P1)', 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q4b_diffusion_reaction.png');

fprintf('Q4b done. n=%d, u_h(0)=%.4f, u_h(1)=%.4f\n', nel, fsol(1), fsol(end));
