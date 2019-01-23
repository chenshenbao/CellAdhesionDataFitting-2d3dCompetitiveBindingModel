clc;
clear;
close all;
%% -----------
global kf2 kr2 kf3 kr3 ac mr ml c m n
tspan = 0:0.01:10;  % s

%  kinetics parameters 
NA = 6.02e23;       % Avogadro constant
kf2 = 1.1e-4;       % 2D on rate,um2/s
kr2 = 7;            % 2D off rate , /s
kf3 = 1.4e5;        % 3D on rate, L/mol.s
kr3 = 0.04;         % 3D off rate, /s
                    % parameters outside contact surface
kf3o = kf3;         % 3D on rate, L/mol.s
kr3o = kr3;         % 3D off rate, /s
j=0.9;
tex = ['kf2=',num2str(kf2),';kr2=',num2str(kr2) ,';kf3=',num2str(kf3), ';kr3=',num2str(kr3)];
xxx = text(1,j,tex);
j=j-0.05;
fname = './fig4B.txt';
z=tspan;
for c0 = 5.5e-8.*[0 1 2 4 8 16]       % soluble molecular density,mol/L

rRBC = 3;             % RBC radius
rTC = 3;              % T cell radius
acl = 4*pi*rTC^2;     % T cell area,um2
acr = 4*pi*rRBC^2;    % RBC area,um2
d = 500e-9;           % contact distance

r0 = 0.35;      % remaining ratio
ks = 0.13/60;   % shedding constant, /s
T = 0;          % shedding time, s
V0 = 1e-3;      % chamber volume, Litre
nc = 1e3;       % cell number
ac = 5;         % contact area,um2
mr0 = 95.8;     % initial receptor density, /um2
ml0 = 122.2;    % initial ligand density,/um2
V = ac*d*1e-9;  % contanct volume,Litre

c = c0+(ml0*acl*nc*(1-r0)*(1-exp(-ks*T)))/(NA*V0); % mol/L

m3o = kf3o*c*acr*mr0./(kf3o*c+kr3o);           % decreasing of psgl-1 number on RBC due to soluble ligand competitive binding 

mr = mr0*kr3o/(kf3o*c+kr3o);                   % psgl-1 density after competitive binding
ml = ml0*(r0+(1-r0)*exp(-ks*T));               % ligand number on T cell after shedding

m = ceil(min(ac*mr,c*V*NA+eps));               % upper limit of 3D bonds
n = ceil(ac*min(mr,ml));                       % upper limit of 2D bonds

% initial condition
p0 = zeros((m+1)*(n+1),1);
p0(1)=1;                                       % no bonds at initial condition

% use ode45 to solve equations
[t,p]=ode45('dpnm',tspan,p0);

% caculate cell adhesion probability (Pmn,n¡Ù0)
p_adh = 1-sum(p(:,1:m+1),2);
z = [tspan;p_adh'];

plot(t,p_adh,'-.','LineWidth',2);drawnow; hold on;  
tex = ['mr=', num2str(mr), ';ml=',num2str(ml),'m3o=', num2str(m3o),'m=', num2str(m), ';n=',num2str(n),';Pa=', num2str(p_adh(end))];
xxx = text(1,j,tex);
j=j-0.1;
end
xlabel('\bf Contact Duration (s)');ylabel('\bf Cell Adhesion Probability');ylim([0 1]);

% export result
resultID=fopen(fname,'a');
fmt=[repmat('%.2d ',1,size(z,1)),'\n'];
fprintf(resultID,fmt,z);
fclose(resultID);
saveas(gcf,'fig4b.fig');
