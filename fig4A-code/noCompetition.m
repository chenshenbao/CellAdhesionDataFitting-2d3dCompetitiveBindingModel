clear;close all;clc;
% for the case considering 2d binding only 
global ac mr ml kf kr n

% set parameters 
ac = 5; mr = 95.8; ml = 122.2; 
kf = 1.1e-4; 
kr = 7;
n = min(ac*[mr,ml]);
tspan = 0:0.01:10; 

% initial condition: no bonds
p0 = zeros(n+1,1);
p0(1) = 1;

% use ode45 to solve the equations here 
[t,p]=ode45('dpn',tspan,p0);

% calculate cell adhesion probability p_adh 
adh = p(:,2:end);
p_adh = sum(adh,2); %(sum(Pmn,n¡Ù0))

% export result 
z=[t,p_adh];
fname = './D2_c0.txt';
resultID=fopen(fname,'a');
fprintf(resultID,'%.2d %.2d\n',z');
fclose(resultID);

% plot 
plot(t,sum(adh,2),'-.','LineWidth',2);hold on;
legends=['kr_2=',num2str(kr,'%1.2e')];
drawnow
ylim([0 1.1]);
axs = gca;
text(1,0.7*axs.YLim(2),{['Ac = ',num2str(ac),'; mr = ',num2str(mr),'; ml = ',num2str(ml)];['kf_2 = ', num2str(kf,'%1.2e'),'; kr_2 = ', num2str(kr,'%1.2e'),'; n = ', num2str(n)]},'FontSize',12)
xlabel('\bf Contact time (s)');ylabel('\bf Cell adhesion probability');
legend(legends,'location','southeast');
title('\bf2D model - parameter test - kr')
saveas(gcf,'D2_c0','fig');