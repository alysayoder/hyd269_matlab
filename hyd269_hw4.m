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
r=data(15); % radial distance from pumping well to observation well
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
Gdistal = ((T./4.064).*406.4)./Limag; 
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
Gintra = ((T./4.064).*406.4)./dx; 
% Gintra is actually just K because the x sec area is equal to the distance
% between nodes 
G = cons*Gintra*-1; 
for n = 2:nnode-1
   G(n,n) = (G(n+1,n)+G(n-1,n)+G(n,n+1)+G(n,n-1)); % assigns G values to each node as the sum of the connections around it
end
% now need to overwrite perimeter nodes affected by gen head BC with G11'=G11-G10
Gdiag = zeros(nnode,nnode);
for i = 1:length(perim)
Gdiag(perim(i),perim(i)) = Gintra-Gdistal; 
end
Gdiag = Gdiag*eye([nnode nnode]);
G = G + Gdiag; 
%% create D matrix
Dval = ((dx.^2).*4.064).*(.00049);
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
numDD = ones(nnode,length(time))*Ho; % this is the mesh populated with initial head everywhere
numDD = (numDD-H).*-1; % subtracting head after pumping from initial head
%% Run at Steady State
%% Compute Mass Balance
% %% Theis analytical solution Jacob
% %r = 300; % distance from node of interest to well
% Q = prate; % m^3 d^-1
% ddVal = zeros(length(time),1);
% for b = 1:length(time)
% dd = (prate./(4.*pi.*T)).*log((2.25.*T.*time(b))./((r.^2).*S)); 
% ddVal(b) = dd;
% end
%% Theis test w/Integral
% well is placed in the upper left corner (1,10)
% verification node is in the center of the grid (5,5)
% triangle sides are 4 and 5, used pythag to get r
%r = 6.403124237432849;
%r = 300; % distance from node of interest to well
Q = prate; % m^3 d^-1
uVal = zeros(1,tmax);
for t = 1:length(time)       % create uVal, an array containing u from 1-100, steps of .01
uVal(t) = (r.^2.*S)./(4.*T.*time(t));
end
%this is the well function evaluated at the u values
fun = @(x) (exp(-x))./x; 
%Populating an array with well function values for each time step. 
WuVal = zeros(1,tmax);
for v = 1:numel(uVal)      
Wu = integral(fun,uVal(v),inf, 'arrayValued', true);
WuVal(v) = Wu; %these make sense w/appendix 1 vals in Fetter
end
%Calculating drawdown at each time step.
ddValI = zeros(1,tmax);
for b = 1:numel(WuVal)
dd = (Q./(4.*pi.*T)).*WuVal(b);
ddValI(b) = dd;
end
% %% Theis Factorial
% %r = 300; % distance from node of interest to well
% Q = prate; % m^3 d^-1
% uVal = zeros(1,tmax);
% for t = 1:length(time)       % create uVal, an array containing u from 1-100, steps of .01
% uVal(t) = (r.^2.*S)./(4.*T.*(time(t)));
% end
% %Populating an array with well function values for each time step. 
% WuVal = zeros(1,tmax);
% for v = 1:numel(uVal)      
% WuVal(v) = -0.5772-log(uVal(v))+uVal(v)-((uVal(v).^2)./(2.*factorial(2)))+((uVal(v).^3)./(3.*factorial(3)))-((uVal(v).^4)./(4.*factorial(4)));
% %these should make sense w/appendix 1 vals in Fetter
% end
% %Calculating drawdown at each time step.
% ddValF = zeros(1,tmax);
% for b = 1:numel(WuVal)
% ddValF(b) =(Q./(4.*pi.*T)).*WuVal(b);
% end
% %% Theis comparison
% loglog(time(4:end),ddVal(4:end),time(4:end),ddValI(4:end),time(4:end),ddValF(4:end));
% set(gca,'XScale','log','YDir','reverse');
% legend('Jacob','Integral','Factorial');
%% verification plot: drawdown vs. time at a node at least two nodes away
% from the well node. Include Theis and numerical model results. 
plot(time(4:end),ddValI(4:end),time(4:end),numDD(230,4:end));
set(gca,'XScale','log','YDir','reverse');%
ylabel('Drawdown (m)');
title('Verification Plot');
xlabel('Time (d)');
legend('Theis','Numerical')
%% surface plot for troubleshooting
figure 
test = mesh;
for n = 1:400
   test(n) = H(n,500);
end
subplot(3,2,1);
surf(test);
axis([0 20 0 20 -5000 1000]);

subplot(3,2,2);
test=mesh;
for n = 1:400
   test(n) = H(n,1000);
end
surf(test);
axis([0 20 0 20 -5000 1000]);

subplot(3,2,3);
test=mesh;
for n = 1:400
   test(n) = H(n,3000);
end
surf(test);
axis([0 20 0 20 -5000 1000]);

subplot(3,2,4);
test=mesh;
for n = 1:400
   test(n) = H(n,6000);
end
surf(test);
axis([0 20 0 20 -5000 1000]);

subplot(3,2,5);
test=mesh;
for n = 1:400
   test(n) = H(n,8000);
end
surf(test);
axis([0 20 0 20 -5000 1000]);

subplot(3,2,6);
test=mesh;
for n = 1:400
   test(n) = H(n,10000);
end
surf(test);
axis([0 20 0 20 -5000 1000]);