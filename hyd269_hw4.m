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

% 24 is the number of connections between nodes
hOld = zeros(1,24);
hNew = {};
P = cell(24);
G = cell(24);
B = {};
X = {};
Y = {};
S = 0.002; 
T = 0.02;
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




