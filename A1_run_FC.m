%%  Flower Code - A two-domain AAA-least squares algorithm for computting the electric (e-)field around a flower
% Includes uniform e-field and arthropods (interior spiders and exterior bees)
% Requires the chebfun package: download from https://www.chebfun.org/
%
% Code:
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off','all'); % close and clear any previous data and turn off warnings

%% Parameter Setup (User Input)
% Switches
VEswit = 1; % plot switch: 0 = plot e-potential, 1 = plot e-field magnitude
COswit = 0; % contour (V) switch: 0 = no contour lines, 1 = plot contour lines
errswit = 0; % error switch: 0 = no error plots, 1 = error with exact (only valid for circular flower bdry, or ellipsodal bdry and no arthropods)

% Flower shape parameters
cswit = 0; % flower switch: 0 = smooth flower, 1 = pointy flower, 
Pts = 500; % no. of points on the flower boundary
npet = 5; % no. of petals
L = 0.01; % petal length [metres]
ang = pi/npet; % angle of flower [radians]
A=1; B=0.0; % POINTY flower: A=1, larger B => more pointy (B=0 for circle)
C = 0.7; % SMOOTH flower: C = indent size (C=0 for circle)

% Background e-field parameters
Einf=0; % Einf = far-field e-strength (set Einf=0 for no uniform stream) [V/m]
eps0 = 8e-12; % permittivity of free space [C/(Vm)]
sus1 = 0; eps1 = 1+sus1; % flower (region 1): e-suscept., rel permittivity [no units]
sus2 = 0; eps2 = 1+sus2; % air (region 2): e-suscept., rel permittivity [no units]

% Arthropod parameters: lambda = charge [C/m] (0 = no charge, >0 = source, <0 = sink), (xp,yp) = coordinates
lambmat = {}; xpmat = {}; ypmat = {}; % cell arrays for input data (comment/uncomment below as desired).
lambmat{1} = 1; xpmat{1} = -1.25; ypmat{1} = 0; % point 1
%lambmat{2} = 1; xpmat{2} = 1.8; ypmat{2} = 1.0; % point 2
%lambmat{3} = 1; xpmat{3} = -1.5; ypmat{3} = -0.7; % point 3 (add more points as desired).

% Mesh grid parameters
meswit = 1; % mesh switch: 0 = x-y meshgrid, 1 = r-theta meshgrid
M = 600; N = 600; MEmax = 4; % mesh spacing in x|r (M) and y|theta (N), max value

%% AAA-LS and Meshgrid Datasets
% Non-dimensional quantities
z = FC_shape(cswit,Pts,npet,ang,A,B,C); % flower shape (petal length 1)
[tileps,efparam,Zpmat,Ajmat] = FC_nondim(L,Einf,lambmat,xpmat,ypmat,z,eps0,eps1,eps2); % nondim params

% Datasets
[xdt, ydt, zdt, xspc, yspc] = FC_mesh(M,N,MEmax,meswit); % complex mesh pts
tic; Va = FC_aaa(z,tileps,efparam,Zpmat,Ajmat); timea = toc; % two-domain AAA-least squares e-potential function handle; timea = time taken for AAA-LS to run
[AAAV, AAAdt, maxAAA, EvAAA] = FC_data(zdt,Va,VEswit,meswit,xspc,yspc); % AAA-LS e-potential ftn handle applied to the mesh

if errswit==0 % ie no comparison with exact solution
    Vdata = {AAAV}; VEdata = {AAAdt}; maxVE = {maxAAA}; Evec = {EvAAA}; % AAA-LS data only
else % comparison with exact solution 
    tic; Ve = FC_exa(z,tileps,efparam,Zpmat,Ajmat); timee = toc; % exact e-potential function handle
    [EXAV, EXAdt, maxEXA, EvEXA] = FC_data(zdt,Ve,VEswit,meswit,xspc,yspc); % exact e-potential ftn handle applied to the mesh
    Vdata = {AAAV,EXAV}; VEdata = {AAAdt, EXAdt}; maxVE = {maxAAA, maxEXA}; Evec = {EvAAA, EvEXA}; % AAA-LS and exact data
end

%% Relative Error with exact solution (only if errswit = 1)
if errswit == 0 % ie no comparison with exact solution
    errdata = 0; maxerr=0; % set these variables to 0
    fprintf("AAA-LS: runtime = " +num2str(timea)+" seconds.\n") % output runtime of two domain AAA-LS
else % comparison with exact solution 
    etoly = 1e-4; % when to switch from relative to absolute error
    [ECe, maxCe, minCe, meanCe] = FC_err(EXAdt,AAAdt,etoly); % comparing exact solutions to AAA-LS
    errdata = {ECe}; maxerr = {maxCe}; % error data across the meshgrid; maximum value of the error
    fprintf("AAA-LS: runtime = " +num2str(timea)+" seconds. Errors: max = "+num2str(maxCe)+", min = "+num2str(minCe)+", mean = "+num2str(meanCe)+".\n") % output runtime and error information
end

%% Plots - contour plot of the e-potential/e-field magnitude across the meshgrid
FC_plot(Vdata,VEdata,maxVE,Evec,errdata,maxerr,VEswit,COswit,errswit,z,efparam,Zpmat,Ajmat,xdt,ydt,1); % plotting function
