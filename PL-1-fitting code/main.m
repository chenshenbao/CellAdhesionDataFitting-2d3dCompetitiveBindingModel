
close all;
clear;
clc;

%% 
global NA data tseq k0 SE ac c  ml0 mr0 V acr ml mr m n kf2 kr2 kf3 kr3 kappa

NA = 6.02e23; % Avogadro constant
%% input experimental data, 
% time series
tt = {
    [0.25	0.5	0.75	1	2	4	6	8	10]
    [0.25	0.5	0.75	1	2	4	6	8	10]
    [0.25	0.5	0.75	1	2	4	6	8	10]
    };
% adhesion probability data
dd = {
   [0.05	0.12	0.25	0.19	0.24	0.29	0.22	0.35	0.26]
   [0.09	0.13	0.17	0.22	0.20	0.27	0.25	0.28	0.25]
   [0.11	0.13	0.18	0.20	0.19	0.19	0.14	0.21	0.15]
    };
% SE of adhesion data
 ss = {
     [0.02	0.02	0.03	0.02	0.03	0.03	0.05	0.03	0.04]
     [0.01	0.02	0.02	0.02	0.02	0.04	0.04	0.04	0.02]
     [0.02	0.02	0.04	0.03	0.03	0.04	0.02	0.01	0.03]};
% soluble PL1 density, Mol/L
 cc = [0;4e-9;4e-8];

% molecule density,[mr,ml], molecules/um2
mm = [64 50;
      64 50;
      64 50];
% contact surface area,um2
aa = [5;5;5];

% Cell parameters
rRBC = 3;
rTC = 3;
acl = 4*pi*rTC^2;     % T cell surface area
acr = 4*pi*rRBC^2;    % RBC surface area
d = 500e-9;           % estimated contact distance
%% 
% set initial and searching range of reaction parameters 
%     kf2 -- kr2 --- kf3 -- kr3 %
k1 = [1   ,   1,     1 ,      1   ]; % initial value for para searching
lb = [1e-8,   1e-1,  1e1,     1e-6]; % lower limits for para searching
ub = [1e-2,   1e2,   1e6,     1e2 ]; % upper limits for para searching

fout = 'PL1-result.txt'; % result file name
        fid = fopen(fout,'a');
        fprintf(fid,'mol c mr ml m n kf2_0 kr2_0 kf3_0 kr3_0 kf2 kr2 kf3 kr3 kappa^2\n');
        fclose(fid);
%% 
for jj = 1:size(tt,1)
tseq0 = tt(jj);
tseq = tseq0{1}';
data0 = dd(jj);
data = data0{1}';
SE0 =   ss(jj);
SE = SE0{1}';
ac = aa(jj);
c = cc(jj);           % soluble antibody initial concentration (mol/L)
mr0 = mm(jj,1);
ml0 = mm(jj,2);
V = ac*d*1e-9;        % space volumn between contact area£¬in unit of Litre
figure(jj);set(gcf,'outerposition',get(0,'screensize'));
spi=1;


for i = -5:-4         % set different initial value of kf2
    for j = -1:0      % set different initial value of kf3
        for x = 4
            for y = -2
                k0 = k1.*[10^(i), 10^(j), 10^(x), 10^(y)];
                figure(jj);
                subplot(2,2,spi);
                
                % call  DeltaP.m to calculate difference between theretical  probability value and experimental value,use lsqnonlin() to estimate the best parameters.
                options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt' ,'Display','iter');
                k = lsqnonlin('DeltaP',k0,lb,ub);
                
                fid = fopen(fout,'a');
                fprintf(fid,'PL1 %.1e %.1e %.1e %.1e %f %f %f %f %f %f %f %f %f %f\n',c,mr,ml,m,n , k0(1),k0(2),k0(3),k0(4),kf2,kr2,kf3,kr3,kappa);
                fclose(fid);
                    
                    spi = spi+1;
            end
        end
    end
end
% show search range on figure
tex = { ['\bf- kf2 -'],['\bf- kr2 -'],['\bf- kf3 -'],['\bf- kr3 -']
        [num2str(ub(1),'%.0e') ],[num2str(ub(2),'%.0e') ],[num2str(ub(3),'%.0e') ],[num2str(ub(4),'%.0e') ]
        [num2str(lb(1),'%.0e') ],[num2str(lb(2),'%.0e') ],[num2str(lb(3),'%.0e') ],[num2str(lb(4),'%.0e') ]
       };
        
[texN,texM]=size(tex);
tmp = get(gca);                                             
x = tmp.XLim;  dx = x(2)-x(1); y = tmp.YLim;  dy = y(2)-y(1);   
texX = x(1) + dx*linspace(0.68,0.95,texM);texY = y(1) + dy*linspace(0.7,0.5,texN);
[posx,posy] = meshgrid(texX, texY);      % set where to put notes 
text(posx(:),posy(:),tex,'Fontsize',9);  % show notes
% show title
xxx = text(-12,4.1,['\bfPL1-C = ',num2str(c),'M']);
saveas(gcf,['\bfPL1-C = ',num2str(c),'M.fig']);
end
