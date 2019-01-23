function dp = dpnm(t,p)
% function dpnm defines the differential equations of 2D and 3D competitive binding:
% dp(n,m)/dt = kf2*ml*[ac*mr-m-(n-1)]*p(n-1,m) + kr2*(n+1).*p(n+1,m) + ...
% kf3*c*[ac*mr-(m-1)-n]*p(n,m-1) + kr3*(m+1).*p(n,m+1) - ...
% [kf2*ml*(ac*mr-m-n)+kr2*n + kf3*c*(ac*mr-m-n)+kr3*m]*p(n,m)
% p(m,n): probability of existing m 3D bonds and n 2D bonds simutaneously, 
% n: number of 2D bonds, 
% m: number of 3D bonds,
% t: reaction time (contact time)
% kf,kr: on rates and off rates

% Chen Shenbao,2017-04, chenshenbao@imech.ac.cn

%% neighbor probability transmission model£º
%    -0--------n-1---------n---------n+1----------------> 2D bond number
%  | -----|-----------|---------|-----------|------
% m-1     |           | p(n,m-1)|           | 
%  | -----|-----------|---------|-----------|------
%  m      |  p(n-1,m) | p(n,m)  |  p(n+1,m) |
%  | -----|-----------|---------|-----------|------
% m+1     |           | p(n,m+1)|           |
% \|/-----|-----------|---------|-----------|------
%  3D bond number

%% boundary conditions
% equations on and outside the boundaries need to be modified. We can find
% the boundaries based on the following 5 possible cases and adjust the
% equations by coefficient matrices as following.
%--------------------------------------------------------------------------
%--------------------case 1------------------------------------------------
%  R: number of membrane receptor            \.R.R.R.R.R.R./   Cell 1
% mL: number of membrane ligand                  sL  sL         space
% sL: number of soluble ligand                .mL.mL.mL.mL.    Cell 2
%                                            /             \
% case 1, R >= (sL+L)      
% or      R >= (m+n)
%    -0--------1---------2---------3--------n-------> 2D bond number
%  | -----|---------|---------|---------|------
%  0  p   |    p    |    p    |    p    |  p    
%  | -----|---------|---------|---------|------
%  1  p   |    p    |    p    |    p    |  p
%  | -----|---------|---------|---------|------
%  m  p   |    p    |    p    |    p    |  p
% \|/-----|---------|---------|---------|------
%  3D bond number
%--------------------------------------------------------------------------
%---------------------case 2-----------------------------------------------
% R: number of receptor            \.R.R.R.R.R.R./
% mL:number of membrane ligand        sL sL sL 
% sL:number of soluble ligand       .mL.mL.mL.mL.
% case 2, m<n<=R<(m+n)             /             \
%
%    -0--------1---------2---------3--------n-------> 2D bond number
%  | -----|---------|---------|---------|------
%0 |  p   |    p    |    p    |    p    |  p    
%  | -----|---------|---------|---------|------
%1 |  p   |    p    |    p    |    p    |  /
%  | -----|---------|---------|---------|------
%m |  p   |    p    |    p    |    /    |  /
% \|/-----|---------|---------|---------|------
%  3D bond number
%--------------------------------------------------------------------------
%--------------------case 3------------------------------------------------
% R: number of receptor            \.R.R.R.R.R.R./
% mL:number of membrane ligand        sL sL sL 
% sL:number of soluble ligand       .mL.mL.mL.mL.
% case 1, max(sL,L)<=R<(sL+L),n<m  /             \
%
%    -0------1------2-------3------n------> 2D bond number
%  | -----|------|------|------|-----|
%0 |  p   |  p   |   p  |   p  |  /  |  
%  | -----|------|------|------|-----|
%1 |  p   |  p   |   p  |   /  |  /  |
%  | -----|------|------|------|-----|
%m |  p   |  p   |   /  |   /  |  /  |
% \|/-----|------|------|------|-----|
%  3D bond number
%--------------------------------------------------------------------------
%------------------case 4---------------------------------------------------
%
%   --0-------1-------2--------> 2D bond number
%  |------|-------|-------|
%0 |  p   |    p  |   p   |    
%  |------|-------|-------|
%1 |  p   |    p  |   p   |   
%  |------|-------|-------|
%2 |  p   |    p  |   p   |   
%  |------|-------|-------|
%3 |  p   |    p  |   /   |      
%  |------|-------|-------|
%m |  p   |    /  |   /   |        
%  |------|-------|-------|
%  m£¬3D bond number
%--------------------------------------------------------------------------
%------------------case 5---------------------------------------------------
%
%   --0------1-------2--------> 2D bond number
%  |------|------|-------|
%0 |  p   |   p  |   p   |
%  |------|------|-------|
%1 |  p   |   p  |   /   |
%  |------|------|-------|
%2 |  p   |   /  |   /   |
%  |------|------|-------|
% \|/
%  m£¬3D bond number
%% 
%---------------------------------------------------------------------------
global kf2 kr2 kf3 kr3 ac mr ml c m n   % global parameters, from main.m
%--------------------------------------------------------------------------
R = ceil(ac*mr);% receptor number

%% ---------------------------------------------------------------
% define differential equations according to regions
[x, y] = meshgrid(0:n, 0:m);% x y correspond to bond number,p(i)=p(i[x,y])
iL = ~(x==0|x+y> R);    iU = ~(y==0|x+y> R); 
iR = ~(x==n|x+y>=R);    iD = ~(y==m|x+y>=R); % coefficient matrix 

i = reshape(1:(m+1)*(n+1),m+1,n+1); % matrix of equantion index
il = i-(m+1); ir = i+(m+1); iu = i-1; id = i+1; % index matrix of leftside/rightside/upside/downside of x
[il(:,1), ir(:,end), iu(1,:), id(end,:)] = deal(1); % deal with the outside region

kL = kf2*ml*(ac*mr-y-(x-1)).*iL; % coefficient matrix 
kR = kr2*(x+1).*iR;
kU = kf3*c*(ac*mr-(y-1)-x).*iU;
kD = kr3*(y+1).*iD;
kI = -(   kf2*ml*(ac*mr-y-x) + kr2*x + kf3*c*(ac*mr-y-x) + kr3*y   );

dp = kL.*p(il) + kR.*p(ir) + kU.*p(iu) + kD.*p(id) + kI.*p(i); % equantion in matrix
dp = dp(:);