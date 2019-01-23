clc;
clear;
close all;
%% -----------
global kf2 kr2 kf3 kr3 ac mr ml c m n
tspan = 0:0.01:10;  % unit, s

%  kinetics parameters 
NA = 6.02e23;       % Avogadro constant
kf2 = 1.1e-4;       % 2D on rate,um2/s
kr2 = 7;            % 2D off rate , /s
kf3 = 1.4e5;        % 3D on rate, L/mol.s
kr3 = 0.04;         % 3D off rate, /s
                    % parameters outside contact surface
kf3o = kf3;         % 3D on rate, L/mol.s
kr3o = kr3;         % 3D off rate, /s
c0 =  5.5e-8;       % 3D molecular density,mol/L

rRBC = 3; % RBC radius
rTC = 3;  % T cell radius
acl = 4*pi*rTC^2;     % T cell area,um2
acr = 4*pi*rRBC^2;    % RBC area ,um2
d = 500e-9;           % contact distance

r0 = 0.35;      % ratio after shedding
ks = 0.13/60;   % shedding constant, s^-1
T = 20*60;      % shedding time, 20min
V0 = 1e-3;      % chamber volume, Litre
nc = 1e3;       % cell number
ac = 5;         % contact area£¬um2
mr0 = 95.8;     % initial receptor density£¬/um2
ml0 = 122.2;    % initial ligand density,/um2
V = ac*d*1e-9;  % volume between contanct area, Litre

c = c0+(ml0*acl*nc*(1-r0)*(1-exp(-ks*T)))/(NA*V0); % mol/L

m3o = kf3o*c*acr*mr0*(1-exp(-(kf3o*c+kr3o)*T))./(kf3o*c+kr3o); % psgl-1 number decreasing by soluble L-selectin competitive binding 
mr = mr0-m3o./acr;                             % psgl-1 density on RBC
ml = ml0*(r0+(1-r0)*exp(-ks*T));               % L-selectin density on T cell after shedding

m = ceil(min(ac*mr,c*V*NA)) ;                  % 3d bonds upper limit
n = ceil(ac*min(mr,ml));                       % 2d bonds upper limit
% initial condition
p0 = zeros((m+1)*(n+1),1);
p0(1)=1;                                       % no bonds at initial condition
% ode45 solving equation
[t,p]=ode45('dpnm',tspan,p0);
% cell adhesion probability (Pmn,n¡Ù0)
p_adh = sum(p(:,m+2:end),2);
plot(t,p_adh,'-.','LineWidth',2);hold on;     

xlabel('\bf Contact Duration (s)');ylabel('\bf Cell Adhesion Probability');ylim([0 1]);
title('23D-entire-shedding');
% notes on figure
tex = {
        [ 'kf2=',  num2str(kf2,'%1.2e') ], [ 'kr2=',  num2str(kr2,'%1.2e') ]
        [ 'kf3=',  num2str(kf3,'%1.2e') ], [ 'kr3=', num2str(kr3,'%1.2e')  ]
        ['c0 = ',num2str(c0,'%1.4e')    ],[ 'c = ',num2str(c,'%1.4e')]
        ['ml = ',num2str(ml,3),' \mum^{-2}'],['mr = ',num2str(mr),'\mum^{-2}']
        [ 'Ac = ',num2str(ac),'\mum^2'],['n = ',num2str(n),';m = ',num2str(m) ]};
[texN,texM]=size(tex);
tmp = get(gca);                                             
x = tmp.XLim;  dx = x(2)-x(1); y = tmp.YLim;  dy = y(2)-y(1);   
if p_adh(end)>=0.5
    texX = x(1) + dx*linspace(0.2,0.4,texM);texY = y(1) + dy*linspace(p_adh(end)*0.75,0.1,texN);
else
    texX = x(1) + dx*linspace(0.05,0.35,texM);texY = y(1) + dy*linspace(0.95,p_adh(end)*1.4,texN);
end
[posx,posy] = meshgrid(texX, texY);  
text(posx(:),posy(:),tex,'Fontsize',10);           

%% 
h1=findobj(gcf,'type','line');
saveas(gcf,'23D_entire_shedding');
x = get(h1,'xdata');
y = get(h1,'ydata');
z = [x',y'];
fname = './23D_entire_shedding.txt';
resultID=fopen(fname,'a');
fprintf(resultID,'%.2d %.2d\n',z');
fclose(resultID);