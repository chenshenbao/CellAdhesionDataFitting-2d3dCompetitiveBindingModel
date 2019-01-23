function dp = dpn(t,p)

%PNM creates the differential equations of 2D  small system.
% probability dynamics model, as follow, 
% dp(n)/dt = ac*mr*ml*kf*p(n-1) - (ac*mr*ml*kf + n*kr)*p(n) + (n+1)*kr*p(n+1)
% p(n): probability of existing n bonds
% n: number of 2D bonds
% t: contact time
%   ----------|------|--------|------|--------|-----|-------|
%         p(1)|  ... | p(i-1) | p(i) | p(i+1) | ... | p(n+1)|
%   ----------|------|--------|------|--------|-----|-------|
% bondnum:  0                    i-1                    n
% Chen Shenbao,2017-04, chenshenbao@imech.ac.cn

global ac mr ml kf kr n
dp = zeros(n+1,1);
x = 1;       dp(x) =  -( ac*mr*ml*kf + (x-1)*kr )*p(x) + x*kr*p(x+1); %n=0
for x = 2:n; dp(x) = ac*mr*ml*kf*p(x-1) - ( ac*mr*ml*kf + (x-1)*kr )*p(x) + x*kr*p(x+1); end
x = n+1;     dp(x) = ac*mr*ml*kf*p(x-1) - ( ac*mr*ml*kf + (x-1)*kr )*p(x); %n=nmax
end
