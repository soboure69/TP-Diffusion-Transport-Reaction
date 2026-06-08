%----------------------------------------------------------------------------
%            
% Solve the ordinary differential equation given as            
%   a u'' + b u' + c u = 1,  0 < x < 1                                                                        
%   u(0) = 0  and  u(1) = 1
% using "nel" linear elements
%
% Variable descriptions                                                      
%   k = element matrix                                             
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
nel=10;                  % number of elements
nnel=2;                 % number of nodes per element
ndof=1;                 % number of dofs per node
nnode=nel+1;                % total number of nodes in system
sdof=nnode*ndof;        % total system dofs  

%-----------------------------------------
%  input data for nodal coordinate values
%-----------------------------------------

h=1/nel;
for i = 1:nnode
    gcoord(i)=(i-1)*h;
end    

%-----------------------------------------------------
%  input data for nodal connectivity for each element
%-----------------------------------------------------

for iel=1:nel           % loop for the total number of elements
    nodes(iel,1)= iel;
    nodes(iel,2)= iel+1;
end

%-----------------------------------------
%  input data for coefficients of the ODE
%-----------------------------------------

acoef=-1;                % coefficient 'a' of the diff eqn
bcoef=0;               % coefficient 'b' of the diff eqn
ccoef=1;                % coefficient 'c' of the diff eqn    

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1)=1;             % first node is constrained
bcval(1)=0;             % whose described value is 0 
bcdof(2)=nnode;             % nnode-th node is constrained
bcval(2)=1;             % whose described value is 1

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % initialization of system force vector
kk=zeros(sdof,sdof);    % initialization of system matrix
index=zeros(nnel*ndof,1);  % initialization of index vector

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

for iel=1:nel           % loop for the total number of elements

nl=nodes(iel,1); nr=nodes(iel,2); % extract nodes for (iel)-th element
xl=gcoord(nl); xr=gcoord(nr);% extract nodal coord values for the element
eleng=xr-xl;            % element length
index=feeldof1(iel,nnel,ndof);% extract system dofs associated with element

k=feode2l(acoef,bcoef,ccoef,eleng); % compute element matrix
f=fef1l(xl,xr);                     % compute element vector
[kk,ff]=feasmbl2(kk,ff,k,f,index);  % assemble element matrices and vectors

end

%-----------------------------
%   apply boundary conditions
%-----------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%----------------------------
%  solve the matrix equation
%----------------------------

fsol=kk\ff;   

%---------------------
% analytical solution  
%---------------------


%for i=1:nnode
%x=gcoord(i);
%esol(i)= ....     %solution exacte sur la grille de calcul
%end

%------------------------------------
% print both exact and fem solutions
%------------------------------------

num=1:1:sdof;
%store=[num' fsol esol']   % si on a une solution exacte esol
store=[num' fsol]
%------------------------------------
% plot both exact and fem solutions
%------------------------------------

%xplot=linspace(0,1);
%esolplot=......    %solution exacte sur grille fine xplot 
%plot(xplot,esolplot,'r',gcoord,esol,'*b')
%legend('exact solution', 'fem solution')
plot(gcoord,fsol,'*b')
legend('fem solution for diffusion transport pb')
%---------------------------------------------------------------

