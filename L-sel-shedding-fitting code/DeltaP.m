function dP =  DeltaP(k) 
% input: k = [kf2 kr2 kf3 kr3 kr3o], output: dP = ss(100*deltaP)

global kf2 kr2 kf3 kr3 kf3o ac mr ml c m n  tseq data k0 SE V ml0 mr0 r0 c0 ks T acr NA kappa

tspan = 0:10/1000:tseq(end);  % time series(s),interval consistent with the experimental contact time

% reaction kinetic parameters 
kf2 = k(1);        % 2D binding rate, um2/s
kr2 = k(2);        % 2D dissociation rate , /s
kf3 = k(3);        % 3D binding rate, L/mol.s
kr3 = k(4);        % 3D dissociation rate, /s
kf3o = kf3;        % parameters outside contact surface
kr3o = kr3;

m3o = kf3o*c*acr*mr0*(1-exp(-(kf3o*c+kr3o)*T))./(kf3o*c+kr3o); % psgl-1 number decreasing by soluble L-selectin competitive binding outside contact region

mr = (acr*mr0-m3o)/acr;                        % density of psgl-1 on rbc
ml = ml0*(r0+(1-r0)*exp(-ks*T));               % density of L-selectin on T cell after shedding 
m = max(ceil(min(ac*mr,c*V*NA)),1);            % upper limit of 3d bonds number, determined by membrane receptor number (ac*mr) and soluble receptor number (c*V*NA)
n = ceil(ac*min(mr,ml));                       % upper limit of 2d bonds number, determined by membrane molecules number on both sides

% set initial condition

p0 = zeros((m+1)*(n+1),1);
p0(1)=1;                                       % no bonds at time zero

% ode45 solving

[t,p]=ode45('dpnm',tspan,p0);

% calculate cell adhesion probability (Pmn,n¡Ù0)

p_adh = sum(p(:,m+2:end),2);

% plot results
errorbar(tseq,data,SE,'ks','MarkerSize',8);hold on;     % bars of experimental data
plot(t,p_adh,'-.','LineWidth',2);hold off;             


xlabel('\bf Contact time (s)');ylabel('\bf Cell Adhesion Probability');ylim([0 1]);
legends = {'Experimental data','ode45 Fitting'};
legend(legends,'fontsize',10);

% export result
y = p_adh(round(1+tseq./0.01));
kappa = sum(((y-data)./SE).^2);
dP = (y-data)./SE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show initial and final parameters on figure
tex = { [      '\bfInitial parameters'   ],[     '\bfFinal parameters'     ]
        ['kf2 = ',num2str(k0(1),'%1.2e') ],[   num2str(kf2,'%1.2e')        ]
        ['kr2 = ',num2str(k0(2),'%1.2e') ],[   num2str(kr2,'%1.2e')        ]
        ['kf3 = ',num2str(k0(3),'%1.2e') ],[   num2str(kf3,'%1.2e')        ]
        ['kr3 = ',num2str(k0(4),'%1.2e') ],[   num2str(kr3,'%1.2e')        ]
        ['\chi^2 = ',num2str(kappa)],[ 'Ac = ',num2str(ac),'\mum^2']
        ['c0 = ',num2str(c0,'%1.2e')     ],[ 'c = ',num2str(c,'%1.2e')     ]
        ['2D bonds, n = ',num2str(n)     ],['3D bonds, m = ',num2str(m)    ]
        ['ml = ',num2str(ml,3),' \mum^{-2}'],['mr = ',num2str(mr),'\mum^{-2}']};
        
[texN,texM]=size(tex);
tmp = get(gca);                                             
x = tmp.XLim;  dx = x(2)-x(1); y = tmp.YLim;  dy = y(2)-y(1);   
if max(data)>=0.5
    texX = x(1) + dx*linspace(0.2,0.4,texM);texY = y(1) + dy*linspace(max(data)*0.75,0.1,texN);
else
    texX = x(1) + dx*linspace(0.05,0.35,texM);texY = y(1) + dy*linspace(0.95,max(data)*1.4,texN);
end
[posx,posy] = meshgrid(texX, texY);  
text(posx(:),posy(:),tex,'Fontsize',10);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow;

end



