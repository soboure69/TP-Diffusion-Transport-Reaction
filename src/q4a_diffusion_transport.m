%% Q4a - Diffusion-Transport : -u'' + u' = 1, u(0)=0, u(1)=1
%  Mise en application DIRECTE des routines fournies (comme Diff_Trans_Reac1D.m)
%  pour les parametres specifiques de Q4a, puis trace de la solution approchee u_h.
%
%  Forme du template fourni :  a u'' + b u' + c u = 1
%  Ici : a=acoef=-1, b=bcoef=1, c=ccoef=0  =>  -u'' + u' = 1

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
bcoef =  1;   % b
ccoef =  0;   % c (= sigma)

%--- conditions aux limites ---
bcdof(1) = 1;     bcval(1) = 0;     % u(0) = 0
bcdof(2) = nnode; bcval(2) = 1;     % u(1) = 1

%--- initialisation ---
ff = zeros(sdof,1);
kk = zeros(sdof,sdof);

%--- assemblage avec les routines fournies ---
for iel = 1:nel
  nl = nodes(iel,1); nr = nodes(iel,2);
  xl = gcoord(nl);   xr = gcoord(nr);
  eleng = xr - xl;
  index = feeldof1(iel,nnel,ndof);          % ddl de l'element
  k = feode2l(acoef,bcoef,ccoef,eleng);     % matrice elementaire
  f = fef1l(xl,xr);                         % vecteur charge (f=1)
  [kk,ff] = feasmbl2(kk,ff,k,f,index);      % assemblage
end

%--- application des conditions de Dirichlet ---
[kk,ff] = feaplyc2(kk,ff,bcdof,bcval);

%--- resolution ---
fsol = kk \ ff;

%--- representation graphique de la solution approchee u_h ---
figure();
plot(gcoord, fsol, 'ro-', 'MarkerSize', 5, 'LineWidth', 1.2);
xlabel('x'); ylabel('u_h(x)');
title(sprintf('Q4a: solution approchee u_h (-u'''' + u'' = 1, n=%d)', nel));
legend('Solution EF u_h (P1)', 'Location', 'NorthWest');
grid on;
print('-dpng', '-r150', '../figures/q4a_diffusion_transport.png');

fprintf('Q4a done. n=%d, u_h(0)=%.4f, u_h(1)=%.4f\n', nel, fsol(1), fsol(end));
