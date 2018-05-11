%% notes
% Objectives:
% -2D, integrated, finite difference
% -solve transient gw flow eqn with an implicit numerical scheme 
% and a direct matrix solver
% -compare with Theis analytical solution 
data = csvread('FD_input.csv',2,0);
nnode=data(1); % get values from csv file
xnode=data(2);
ynode=data(3);
T = data(4); % transmissivity or K because b=1 here
S = data(5); % storativity
dt = data(6); % time interval
dx = data(7); % distance between each node, also x sec are because b=1
wellnode=data(8);
prate=data(9);
Ho = data(10); % define initial head 
tmax=data(11); % max time value
qp=data(12); % flux at perimter cells, these have gen. head BCs
himag=data(13); % head at imaginary nodes for gen head BC
Limag=data(14); % distance from node to imaginary nodes

%% create domain w/numbered nodes
mesh = reshape(1:nnode, [xnode ynode]);
%% add general head boundary conditions
q = zeros(nnode,1); % flux, this is the boundary condition vector
perim = zeros((xnode.*2)+(ynode.*2)-4,1); % perim holds indexes of all perimeter nodes
for i = 1:xnode
    perim(i)=mesh(i,1);
end
for i= 1:xnode
    perim(i+xnode)=mesh(i,ynode);
end
for j = 2:ynode-1
    perim(j+(2.*xnode)-1)=mesh(1,j);
end
for j = 2:ynode-1
    perim(j+(2.*xnode)+ynode-3)=mesh(xnode,j);
end
% take index of perimeters from above and populate q with G for these
% locations
Gdistal = -(T.*dx)./Limag; 
for w = perim(1:end)
q(w,1)=qp-((Gdistal).*himag); 
end
q(wellnode)=prate./dx/dx; % define pumping rate as specified flux BC
%% create G matrix
% create logical nnode x nnode matrix
% this will be multiplied by the G value to create a matrix of
% G values populating connected node slots
cons = zeros(nnode,nnode);
for i=1:nnode
    for j=1:nnode
        if i == j 
            cons(i,j) = 0;
        elseif j == i+1 || j == i-1
            cons(i,j)=1;
        elseif i == j+1 || i == j-1
            cons(i,j)=1;
        elseif j == i+xnode || j == i-xnode
            cons(i,j)=1;
        elseif i == j+ynode || i == j-ynode
            cons(i,j)=1;
        else
            cons(i,j)=0;
        end
    end
end
% populate G matrix with G values at connected nodes and boundary nodes
Gintra = -(T.*dx)./dx; 
% Gintra is actually just K because the x sec area is equal to the distance
% between nodes 
G = zeros(nnode,nnode); % conductance array (sparse)
G = cons*Gintra; 
% now need to overwrite nodes affected by gen head BC with G11'=G11-G10
Gdiag = eye([nnode nnode])*(Gintra-Gdistal);
G = G + Gdiag; 
%% create D matrix
Dval = (dx.^2).*S;
D = eye([nnode nnode])*Dval; % capacitance array (vector)
%% initialize arrays hOld and hNew
hO = ones(nnode,1)*10; % set initial head to 10 m at all nodes
hOld = ones(nnode,1);
hNew = ones(nnode,1);
% hOld is the head at time step n, hNew is nead at n+1
%% Theis analytical solution
% well is placed in the upper left corner (1,10)
% verification node is in the center of the grid (5,5)
% triangle sides are 4 and 5, used pythag to get r
%r = 6.403124237432849;
r = 200; % distance from node of interest to well
Q = 2000; % m^3 d^-1
uVal = zeros(1,tmax.^2);
time=zeros(tmax,1);
for i=1:tmax
time(i) = time+dt.*1.2;
end
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