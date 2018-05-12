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
prate=data(9); %m^3/day
Ho = data(10); % define initial head 
tmax=data(11); % max time value
qp=data(12); % flux at perimiter cells, these have gen. head BCs
himag=data(13); % head at imaginary nodes for gen head BC
Limag=data(14); % distance from node to imaginary nodes
%% create domain w/numbered nodes
mesh = reshape(1:nnode, [xnode ynode]);
%% create time vector used for Theis and Direct Solver
% time=zeros(tmax,1);
% time_factor = dt; 
% for i=1:tmax
% time(1,1)=1;
% time(1+i,1) = time(i,1)+time_factor;
% time_factor = 1.2.*(time(i+1)-time(i)); 
% end
% time=time./86400;
time = dt:dt:tmax;
%% add general head boundary conditions
q = zeros(nnode,1); % flux, this is the boundary condition vector
perim = zeros((xnode.*2)+((ynode.*2)-4),1); % perim holds indexes of all perimeter nodes
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
q(wellnode)=-prate./dx./dx; % define pumping rate as specified flux BC
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
Gintra = -(.02.*dx)./dx; 
% Gintra is actually just K because the x sec area is equal to the distance
% between nodes 
G = cons*Gintra; 
% now need to overwrite nodes affected by gen head BC with G11'=G11-G10
Gdiag = eye([nnode nnode])*(Gintra-Gdistal);
G = G + Gdiag; 
%% create D matrix
Dval = (dx.^2).*S;
D = eye([nnode nnode])*Dval; % capacitance array (vector)
%% initialize array for storing H values for each node at each time step
hInit = ones(nnode,1)*Ho; 
H = zeros(nnode, length(time));
H(:,1) = hInit; % set initial head to 10 m at all nodes for 1st time step
%% Format Ax=b matrix to be solved at each time step
A = (dt*G+D);
for i = 1:length(time)
B = D*H(find(H,nnode,'last')) + q*time(i); % selects last column of H
x = mldivide(A,B);
H(:,i) = x;
end
%% Compute Numerical Drawdown
numDD = ones(nnode,length(time))*Ho;
numDD = numDD-H;
%% Theis analytical solution
r = 100; % distance from node of interest to well
Q = prate; % m^3 d^-1
ddVal = zeros(length(time),1);
for b = 1:length(time)
dd = (prate./(4.*pi.*T)).*log((2.25.*T.*time(b))./((r.^2).*S)); 
ddVal(b) = dd;
end
%% verification plot: drawdown vs. time at a node at least two nodes away
% from the well node. Include Theis and numerical model results. 
plot(time,ddVal,time,numDD(193,:));
set(gca,'XScale','log','YDir','reverse');
ylabel('Drawdown (m)');
title('Verification Plot');
xlabel('Time (d)');
legend('Theis','Numerical')