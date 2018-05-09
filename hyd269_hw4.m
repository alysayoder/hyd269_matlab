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
[a,b,data2] = xlsread('FD_input.csv');
% Compute G and D matrices, then use them to formulate Ax=b and use
% a direct solver to find the solution at each time step. 
% in fully implicit methods, alpha=1
% 24 is the number of connections between nodes
%% manually entering variables, eventually will be file sourced
% on p. 78 W&A, well is discharging constantly at 2000 m^3 d^-1
hOld = zeros(1,24);
hNew = {};
P = cell(24);
G = cell(24);
B = {};
X = {};
Y = {};
S = 0.002; 
T = 0.02; % on p.78 of the book, T=300 m^2 d^-1
XSI = [-.57735,.57735,.57735,-.57735];
ETA = [-.57735,-.57735,.57735,.57735];
NS = zeros(1,4);
NX = zeros(1,4);
NY = zeros(1,4);
node = zeros(1,4);
%% block 1: generate nodal coordinates & initial & boundary conditions
% 16 is the number of nodes and elements, this will be in the data file
nnode = 16;
nelem = 16; % may not need this for FD
%% Theis analytical solution
% well is placed in the upper left corner (1,10)
% verification node is in the center of the grid (5,5)
% triangle sides are 4 and 5, used pythag to get r
r = 6.403124237432849;
Q = 2000; % m^3 d^-1
uVal = zeros(1,100);
for time = 1:100         % create uVal, an array containing u from 1-100
u = (r.^2.*S)./(4.*T.*time);
uVal(time) = u;
end
% this is the well function evaluated at the u values
fun = @(uVal) (exp(-uVal(i)))./uVal(i); 
% useq = zeros(1,100);
% for j = 1:numel(uVal)
%     po = fun(uVal(j));
%     useq(j) = po;
% end
%%
WuVal = zeros(1,100);
for i = 1:numel(uVal)      
Wu = integral(fun,uVal(i),inf, 'arrayValued', true);
WuVal(i) = Wu; %these make sense w/appendix 1 vals in Fetter
end

ddVal = zeros(1,100);
for r = 1:numel(WuVal)
dd = (Q./(4.*pi.*T)).*WuVal(r);
ddVal(r) = dd;
end
%% verification plot: drawdown vs. time at a node at least two nodes away
% from the well node. Include Theis and numerical model results. 






