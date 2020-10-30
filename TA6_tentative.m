
global betta theta r nj nx ny pi gridy pens sr curv grdfac jr lphi

%% calibration 

r = .045; 
rho = .05; 
betta = 1/(1+rho); 
theta = 2; 
lphi = .5; 

nj = 70; 
ny = 2; 
nx = 40; 
jr = 45;

% Grids
% two-state markov process
gridy = [0.5; 1.5];
pi = [.6, .4; .4, .6]; 

% deterministic income
w = 1; 
wprofile = ones(nj,1); 
tau = 0; 
netw = w*wprofile*(1-tau);
iw = ones(nj,1);
iw(jr:nj) = 0;
pens = 0; 

% grid saving 
curv = 3;       % curvature of grid
grdfac = 40;      % scaling factor of saving grid


%% solution

disp('solution of household model');

% grids and decisions rules:
gridx = zeros(nj,ny,nx);
gridsav = zeros(nx,1);
gridass = zeros(nj,ny,nx);
cfun = zeros(nj,ny,nx);
vfun = zeros(nj,ny,nx);
vpfun = zeros(nx,1);
vptrans = zeros(nj,ny,nx);
lfun = zeros(nj,ny,nx); 

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(2:nx)=makegrid(0.0,grdfac,nx-1,curv);
gridsav(1)=0.0;

% final period
for yc=1:ny
% cash-on-hand grid at nj:
inc = iw(nj)*netw(nj)*gridy(yc)+(1-iw(nj))*pens;
% in case of no pension system, assume some minimum cash on hand:
minx=max(inc,sqrt(eps));
maxx=gridsav(nx)*(1.0+r)+inc;
gridx(nj,yc,:)=linspace(minx,maxx,nx);
% Final period Consumption function, asset holdings, value function, including derivative
cfun(nj,yc,:)=gridx(nj,yc,:); 
lfun(nj,yc,:)=0; 
gridass(nj,yc,:)=(gridx(nj,yc,:)-inc)/(1+r);
vfun(nj,yc,:)= U(cfun(nj,yc,:),1); % utility function
vpfun(:) = UC(cfun(nj,yc,:),1);  % marginal utility wrt c
vptrans(nj,yc,:)= vpfun.^(1/(lphi*(1-theta)-1));  
end

% iterate backwards at retirement 
% ---------------------------------------------------------------------------------------------------------------------------------------------
for jc=nj-1:-1:jr
for yc=1:ny    
for xc=2:nx    % no binding borrowing limit x > xmin
vp=zeros(2,1);  
lfun(jc,yc,:) = zeros(nx,1);

for ycc=1:ny
% income tomorrow:
incp1=iw(jc+1)*netw(jc+1)*gridy(ycc)+(1-iw(jc+1))*pens;
% Maximum cash on hand tomorrow:
% in case of zero savings and no pension system assume some
% minimum cash on hand
cah=max(sqrt(eps),incp1+(1.0+r)*gridsav(xc));
% Interpolate derivative of value function
if (cah<gridx(jc+1,ycc,1))
disp('how can this be?')
end
if ( cah>gridx(jc+1,ycc,nx) )
% if out of bounds simply set it to decision at nx:
vptr = vptrans(jc+1,ycc,nx);
else
vptr = interp1(squeeze(gridx(jc+1,ycc,:)),squeeze(vptrans(jc+1,ycc,:)),cah);
end
vp(ycc)=vptr.^(lphi*(1-theta)-1);
end

% Euler equation: RHS
expvp=betta*(1.0+r)*sum(pi(yc,:)*vp(:));
% consumption
cfun(jc,yc,xc)=(expvp).^(1/(lphi*(1-theta)-1));
% endogenous x-grid:
gridx(jc,yc,xc) = gridsav(xc)+cfun(jc,yc,xc);

end  % end x > xmin

% income (wages and pensions) in current period/age:
inc=iw(jc)*netw(jc)*gridy(yc)+(1-iw(jc))*pens;
% decision at minx
% notice: correction required for welfare calculation
% the above is actually slightly inefficient because xmin
% can be explicitly computed, then gridsav would be age and
% state dependent.
minx=max(inc,sqrt(eps));
if (minx<gridx(jc,yc,2))
gridx(jc,yc,1)=minx;
else    % set it to some arbitrary fracion of x(2)
gridx(jc,yc,1)=0.9*gridx(jc,yc,2);
end

% Compute optimal consumption and leisure for minx
cfun(jc,yc,1)=gridx(jc,yc,1);
% assets at all xc:
gridass(jc,yc,:)=(gridx(jc,yc,:)-inc)/(1+r);

% Update vfun and vpfun
vpfun(:)=lphi*(cfun(jc,yc,:).^(lphi*(1-theta)-1)).*((1-lfun(jc,yc,:)).^((1-lphi)*(1-theta)));  % marginal utility wrt c;
vptrans(jc,yc,:)=vpfun(:).^(1/(lphi*(1-theta)-1));

% Calculate value function
for xc=1:nx
v=zeros(2,1);
for ycc=1:ny
% income tomorrow:
incp1=iw(jc+1)*netw(jc+1)*gridy(ycc)+(1-iw(jc+1))*pens;
% cah tomorrow
cah=max(sqrt(eps),incp1+(1.0+r)*gridsav(xc));
% this should never be the case:
if ((cah+0.0001)<gridx(jc+1,ycc,1))
warning('How can this be ?');
end
% linear interpolation:
v(ycc)=func_intp(squeeze(gridx(jc+1,ycc,:)),squeeze(vfun(jc+1,ycc,:)),cah);
end    % end for ycc
% update value function
expv=sum(pi(yc,:)*v(:));
u = (1/(1-theta)).*((cfun(nj,yc,xc).^lphi).*((1-lfun(nj,yc,xc)).^(1-lphi))).^(1-theta);
vfun(jc,yc,xc)= u + betta*expv;
end   % end for xc
end   % end for yc
end   % end for jc

% iterate last period of work
% ----------------------------------------------------------------------------------------------------------------------------------------------------
jc = jr-1;
for yc=1:ny    
for xc=2:nx    % no binding borrowing limit x > xmin
for ycc=1:ny
% income tomorrow:
incp1=iw(jc+1)*netw(jc+1)*gridy(ycc)+(1-iw(jc+1))*pens;
% Maximum cash on hand tomorrow:
% in case of zero savings and no pension system assume some
% minimum cash on hand
cah=max(sqrt(eps),incp1+(1.0+r)*gridsav(xc));
% Interpolate derivative of value function
if ( cah<gridx(jc+1,ycc,1))
disp('how can this be?')
end
if ( cah>gridx(jc+1,ycc,nx) )
% if out of bounds simply set it to decision at nx:
vptr = vptrans(jc+1,ycc,nx);
else
vptr = interp1(squeeze(gridx(jc+1,ycc,:)),squeeze(vptrans(jc+1,ycc,:)),cah);
end
vp(ycc)=vptr.^(lphi(1-theta)-1);
end

% Euler equation: RHS
expvp=betta*sr(jc)*(1.0+r)*sum(pi(yc,:)*vp(:));
% consumption if lfun(jc,yc,:)= zeros(nx,1)
cfun(jc,yc,xc)=(expvp).^(1/(lphi(1-theta)-1));
% find KT multiplier from FOC
mul = (1-lphi).*c(jc,yc,xc).^((1-lphi)*(1-theta)) - expvp.*netw(jc).*gridy(yc);
if mul < 0 
cfun(jc,yc,xc)=((1/lphi).*(((lphi/(1-lphi))*netw(jc).*gridy(yc)).^((1-lphi)*(1-theta))).*expvp).^(-1/theta);
lfun(jc,yc,xc) = 1 - (1/(netw(jc).*gridy(yc)))*((1-lphi)/lphi)*cfun(jc,yc,xc);
end 
% endogenous x-grid:
gridx(jc,yc,xc)=gridsav(xc)+cfun(jc,yc,xc);

end  % end x > xmin

% income (wages and pensions) in current period/age:
inc=epsi(jc)*netw(jc)*gridy(yc)+(1-epsi(jc))*pens;
% decision at minx
% notice: correction required for welfare calculation
% the above is actually slightly inefficient because xmin
% can be explicitly computed, then gridsav would be age and
% state dependent.
minx=max(inc,sqrt(eps));
if (minx<gridx(jc,yc,2))
gridx(jc,yc,1)=minx;
else    % set it to some arbitrary fracion of x(2)
gridx(jc,yc,1)=0.9*gridx(jc,yc,2);
end

% [...]

% if gridx(jc,yc,1) > netw(jc)*gridy(1)
% xadd = [netw(jc)*gridy(1); gridx(jc,yc,1) - 0.001; gridx(jc,yc,:)];
% cfun(jc,yc,xc) = lphi*xadd(jc,yc,xc); 
% cfun(jc,yc,1) = 1 - (1/(netw(jc)*gridy(yc)))*(1-lphi)*xadd(jc,yc,xc);
% end

% Compute optimal consumption and leisure for minx
cfun(jc,yc,1)=gridx(jc,yc,1);
% assets at all xc:
gridass(jc,yc,:)=(gridx(jc,yc,:)-inc)/(1+r);

% Update vfun and vpfun
vpfun(:)=lphi*(cfun(jc,yc,:).^(lphi*(1-theta)-1)).*((1-lfun(jc,yc,:)).^((1-lphi)*(1-theta)));  % marginal utility wrt c;
vptrans(jc,yc,:)=vpfun(:).^(1/(lphi(1-theta)-1));

% Calculate value function
for xc=1:nx
v=zeros(2,1);
for ycc=1:ny
% income tomorrow:
incp1=epsi(jc+1)*netw*gridy(ycc)+(1-epsi(jc+1))*pens;
% cah tomorrow
cah=max(sqrt(eps),incp1+(1.0+r)*gridsav(xc));
% this should never be the case:
if ((cah+0.0001)<gridx(jc+1,ycc,1))
warning('How can this be ?');
end
% linear interpolation:
v(ycc)=func_intp(squeeze(gridx(jc+1,ycc,:)),squeeze(vfun(jc+1,ycc,:)),cah);
end    % end for ycc
% update value function
expv=sum(pi(yc,:)*v(:));
u = (1/(1-theta)).*((cfun(nj,yc,xc).^lphi).*((1-lfun(nj,yc,xc)).^(1-lphi))).^(1-theta);
vfun(jc,yc,xc)= u + betta*sr(jc)*expv;
end   % end for xc
end   % end for yc
% end for jc   


% working age 
% ----------------------------------------------------------------------------------------------------------------------------------------------
for jc = jr-2:-1:1
for yc=1:ny    
for xc=2:nx    % no binding borrowing limit x > xmin
for ycc=1:ny
% income tomorrow:
incp1=iw(jc+1)*netw(jc+1)*gridy(ycc)+(1-iw(jc+1))*pens;
% Maximum cash on hand tomorrow:
% in case of zero savings and no pension system assume some
% minimum cash on hand
cah=max(sqrt(eps),incp1+(1.0+r)*gridsav(xc));
% Interpolate derivative of value function
if ( cah<gridx(jc+1,ycc,1))
disp('how can this be?')
end
if ( cah>gridx(jc+1,ycc,nx) )
% if out of bounds simply set it to decision at nx:
vptr = vptrans(jc+1,ycc,nx);
else
vptr = interp1(squeeze(gridx(jc+1,ycc,:)),squeeze(vptrans(jc+1,ycc,:)),cah);
end
vp(ycc)=vptr.^(lphi(1-theta)-1);
end

% Euler equation: RHS
expvp=betta*sr(jc)*(1.0+r)*sum(pi(yc,:)*vp(:));
% consumption if lfun(jc,yc,:)= zeros(nx,1)
cfun(jc,yc,xc)=(expvp).^(1/(lphi(1-theta)-1));
% find KT multiplier from FOC
mul = (1-lphi).*c(jc,yc,xc).^((1-lphi)*(1-theta)) - expvp.*netw(jc).*gridy(yc);
if mul < 0 
cfun(jc,yc,xc)=((1/lphi).*(((lphi/(1-lphi))*netw(jc).*gridy(yc)).^((1-lphi)*(1-theta))).*expvp).^(-1/theta);
lfun(jc,yc,xc) = 1 - (1/(netw(jc).*gridy(yc)))*((1-lphi)/lphi)*cfun(jc,yc,xc);
end 
% endogenous x-grid:
gridx(jc,yc,xc)=gridsav(xc)+cfun(jc,yc,xc);

end  % end x > xmin

% income (wages and pensions) in current period/age:
inc=epsi(jc)*netw(jc)*gridy(yc)+(1-epsi(jc))*pens;
% decision at minx
% notice: correction required for welfare calculation
% the above is actually slightly inefficient because xmin
% can be explicitly computed, then gridsav would be age and
% state dependent.
minx=max(inc,sqrt(eps));
if (minx<gridx(jc,yc,2))
gridx(jc,yc,1)=minx;
else    % set it to some arbitrary fracion of x(2)
gridx(jc,yc,1)=0.9*gridx(jc,yc,2);
end

% [...]

% if gridx(jc,yc,1) > netw(jc)*gridy(1)
% xadd = [netw(jc)*gridy(1); gridx(jc,yc,1) - 0.001; gridx(jc,yc,:)];
% cfun(jc,yc,xc) = lphi*xadd(jc,yc,xc); 
% cfun(jc,yc,1) = 1 - (1/(netw(jc)*gridy(yc)))*(1-lphi)*xadd(jc,yc,xc);
% end

% Compute optimal consumption and leisure for minx
cfun(jc,yc,1)=gridx(jc,yc,1);
% assets at all xc:
gridass(jc,yc,:)=(gridx(jc,yc,:)-inc)/(1+r);

% Update vfun and vpfun
vpfun(:)=lphi*(cfun(jc,yc,:).^(lphi*(1-theta)-1)).*((1-lfun(jc,yc,:)).^((1-lphi)*(1-theta)));  % marginal utility wrt c;
vptrans(jc,yc,:)=vpfun(:).^(1/(lphi(1-theta)-1));

% Calculate value function
for xc=1:nx
v=zeros(2,1);
for ycc=1:ny
% income tomorrow:
incp1=epsi(jc+1)*netw*gridy(ycc)+(1-epsi(jc+1))*pens;
% cah tomorrow
cah=max(sqrt(eps),incp1+(1.0+r)*gridsav(xc));
% this should never be the case:
if ((cah+0.0001)<gridx(jc+1,ycc,1))
warning('How can this be ?');
end
% linear interpolation:
v(ycc)=func_intp(squeeze(gridx(jc+1,ycc,:)),squeeze(vfun(jc+1,ycc,:)),cah);
end    % end for ycc
% update value function
expv=sum(pi(yc,:)*v(:));
vfun(jc,yc,xc)=U(cfun(jc,yc,xc))+betta*sr(jc)*expv;
end   % end for xc
end   % end for yc
end   % end for jc

%% Functions

function y = U(c,l)

global theta lphi

y = (1/(1-theta)).*(c.^lphi.*(1 - l).^(1-lphi)).^(1-theta);

end

function y = UC(c,l)

global theta lphi

y = lphi*(c.^(lphi*(1-theta)-1)).*((1-l).^((1-lphi)*(1-theta))); 

end

function fv = func_intp(x,func,xp)


n = length(x);
if ( xp>x(n) )
% fv = func(n);
fv=func_extrapol(x(n-1),x(n),func(n-1),func(n),xp);
elseif (xp<x(1))
% fv = func(1);
fv=func_extrapol(x(1),x(2),func(1),func(2),xp);
else
fv = interp1(x,func,xp);
end

end


function y=func_extrapol(x1,x2,y1,y2,x)

% simple linear extrapolation

m = (y2-y1)/(x2-x1);
y = y1 + m*(x-x1);

end
    
    function grd = makegrid(x1,x2,n,c)
% makes curved grid according to curvature parameter c
scale=x2-x1;
grd(1)=x1;
grd(n)=x2;
for i=2:n-1
    grd(i)=x1+scale*((i-1.0)/(n-1.0))^c;
end
end
