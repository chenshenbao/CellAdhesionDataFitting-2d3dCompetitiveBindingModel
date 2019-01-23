clear;close all;clc;
% for the case that considering ligands concentration is c = 5.5e-8 mol/L outside contact area
global ac mr ml kf kr n

% set parameters 
ac = 5;             % contact area,um2
mr0 = 95.8;         % initial receptor density, um^-2
ml0 = 122.2;        % initial ligand density,um^-2
rRBC = 3;           % RBC radius,um
acr = 4*pi*rRBC^2;  % RBC area,um2
kf3 = 1.5e5;        % 3D on rate, L mol^-1 s^-1
kr3 = 4e-2;         % 3D off rate, s^-1
kf3o = 1.4e5;       % 3D off rate, L mol^-1 s^-1
kr3o = 4e-2;        % 3D off rate, s^-1
c = 5.5e-8;         % mol/L
m3o = kf3o*c*acr*mr0./(kf3o*c+kr3o); % psgl-1 number decreasing by soluble L-selectin competitive binding outside contact region
mr = mr0-m3o./acr;                   % psgl-1 density on RBC
ml = ml0;                            % L-selectin density on T cell 


kf = 1.1e-4;        % on rate
kr = 7;             % off rate
n = ceil(min(ac*[mr,ml]));
tspan = 0:0.01:10; % 10s
% initial condition: no bonds
p0 = zeros(n+1,1);
p0(1) = 1;

% use ode45 to solve the equations here 
[t,p] = ode45('dpn',tspan,p0);

% calculate cell adhesion probability p_adh
adh = p(:,2:end);
p_adh = sum(adh,2);

% export result 
z = [t,p_adh];
fname = './D23_c1.txt';
resultID=fopen(fname,'a');
fprintf(resultID,'%.2d %.2d\n',z');
fclose(resultID);

% plot 
plot(t,sum(adh,2),'-.','LineWidth',2);hold on;
legends = ['kr_2=',num2str(kr,'%1.2e')];
drawnow
ylim([0 1.1]);
axs = gca;
text(1,0.45*axs.YLim(2),['Ac = ',num2str(ac),'; mr = ',num2str(mr),'; ml = ',num2str(ml),'; kf_2 = ', num2str(kf,'%1.2e'),';'],'FontSize',12)
xlabel('\bf Contact time (s)');ylabel('\bf Cell adhesion probability');
legend(legends,'location','southeast');
title('\bf2D model - parameter test - kr')
saveas(gcf,'D23_c1','fig');