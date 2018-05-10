%% notes
% Objectives:
% -2D, integrated, finite difference
% -solve transient gw flow eqn with an implicit numerical scheme 
% and a direct matrix solver
% -compare with Theis analytical solution 

% The program should be designed such that node numbers, 
% aquifer properties,node-to-node  connections  (including 
% cross-sectional  areas  and  distances)sources  and  sinks,  
% and  any other boundary condition information are read in 
% from a data file.
data = csvread('FD_input.csv',1,0);
%[a,b,data2] = xlsread('FD_input.csv');
% Compute G and D matrices, then use them to formulate Ax=b and use
% a direct solver to find the solution at each time step. 
% in fully implicit methods, alpha=1
nnode=data(1);
xnode=data(2);
ynode=data(3);
tmax=data(12);
% get values from csv file
%% manually entering variables, eventually will be file sourced
% on p. 78 W&A, well is discharging constantly at 2000 m^3 d^-1
hOld = zeros(xnode,ynode);
hNew = zeros(xnode,ynode);
R = zeros(xnode,ynode);
DD = zeros(xnode,ynode);
D = ones(xnode,ynode);
G = ones(xnode,ynode);
S = data(4); 
T = data(5); % on p.78 of the book, T=300 m^2 d^-1
dt = data(6); 
dx = data(7);
Ho = data(11); % define initial head 
tol=0.001; % may not need this
% hOld is the head at time step n, hNew is nead at n+1
%% initialize arrays
for i=1:23
    for j=1:23
        hNew(i,j)=Ho; % set hNew and hOld equal to 10 in all cells
        hOld(i,j)=Ho;
        R(i,j)=0; % set R matrix to 0 in all cells
    end
end
%% define pumping rate as recharge to cell 
R(data(8),data(9))=data(10)./dx/dx;
ti=0:dt:tmax;
%% dig in: create G and D matrices
% assume thickness (b) = 1
% using gwmodfundamentals handout equation 53
% for i =1:23
%     for j=1:23
%         G(i,j) = -T; % from G(i,j) = -T/b * delta x /delta y
%     end
% end
% 
% for i = 1:23
%     for j = 1:23
%         D(i,j) = S.*dx.^2;
%     end
% end

% modified from W&A p. 82 figure 4.6
alpha = 1;
for it=ti(1:end)
    for i=2:22
    for j=2:22
h1 = (hOld(i,j+1)+hOld(i,j-1)+hOld(i+1,j)+hOld(i-1,j))/4;
h2 = (hNew(i,j+1)+hNew(i,j-1)+hNew(i+1,j)+hNew(i-1,j))/4;
f1 = ((dx.^2).*S)/(4.*T.*dt);
f2 = 1./(f1+alpha);
G(i,j) = f2;
D(i,j) = (f1.*hOld(i,j))+(1-alpha).*(h1-hOld(i,j))+(alpha.*h2)+(R(i,j).*dx.^2)/(4.*T);
    end
end

% add boundary conditions, these are no flows on all perimeters
for i = 2:22
    hNew(i,1)=hNew(i,3);
    hNew(i,23)=hNew(i,21);
end
for j = 2:22
    hNew(1,j)=hNew(3,j);
    hNew(23,j)=hNew(21,j);
end
it=it+1;
end
% solve for hNew
plz = G\D; % matrix G currently has no inverse, can't solve yet
%% Theis analytical solution
% well is placed in the upper left corner (1,10)
% verification node is in the center of the grid (5,5)
% triangle sides are 4 and 5, used pythag to get r
%r = 6.403124237432849;
r = 200; % distance from node of interest to well
Q = 2000; % m^3 d^-1
uVal = zeros(1,tmax.^2);
time = 0:dt:tmax;
for t = 1:length(time)       % create uVal, an array containing u from 1-100, steps of .01
u = (r.^2.*S)./(4.*T.*time(t));
uVal(t) = u;
end
%% this is the well function evaluated at the u values
fun = @(uVal) (exp(-uVal(v)))./uVal(v); 
%% Populating an array with well function values for each time step. 
WuVal = zeros(1,tmax.^2);
for v = 1:numel(uVal)      
Wu = integral(fun,uVal(v),inf, 'arrayValued', true);
WuVal(v) = Wu; %these make sense w/appendix 1 vals in Fetter
end
%% Calculating drawdown at each time step.
ddVal = zeros(1,tmax.^2);
for b = 1:numel(WuVal)
dd = (Q./(4.*pi.*T)).*WuVal(b);
ddVal(b) = dd;
end
%% verification plot: drawdown vs. time at a node at least two nodes away
% from the well node. Include Theis and numerical model results. 
theis = timeseries(ddVal, time);
plot(theis);
set(gca,'XScale','log','YDir','reverse','YScale','log');
ylabel('Drawdown (m)');
title('Verification Plot');
xlabel('Time (d)');


