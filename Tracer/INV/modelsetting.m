function [param,iSol]=modelsetting(nmeas,obsindex,obs_path,param)

% flag for running INIT or INV mode
param.init = 0; 

%% Setting up folders and MPI for each delta
param.r = [param.rank;param.rank;param.rank];
nfolder = max(param.r);
param.nodes = param.rank; 
param.ppn = 24;
% processor per simulation, 1000 grid points per proessor, 29806 / 24 ~ 1241
param.pps = 24; % #processor per simulation
creating_directories(); 
writing_processors_name(param,nfolder); 


%%%%%%% Setup parameters  %%%%%%%%

% Grid information

param.n = nmeas; %80;
load coord_htr2.mat; 
param.x = x_htr2; param.y = y_htr2; param.z = z_htr2; 
param.G0= readMESH('./tmp/state_update/MESH');
param.m = length(param.x); 
param.M = 29806; % the size of the full domain

% state variable noise P
param.x_std = [100000;0.5;2]; % pressure;logk

% dynamic model noise Q
param.dx_std =[100;0.01;0.1]; % [pressure; logk]

% R of the scaled observations [pressure, logk]
%param.obsstd = [0.0001;0.01]; 

param.obsindex = obsindex;
param.nmeas = max(param.obsindex{end});
if (param.n-param.nmeas)~=0 ; disp('Measurement inconsistency'); keyboard; end;

% kernel
param.L = [100;100;10];
% NOTE: h should include correlation lengths 
param.kernel = @(h) exp(-h.^2); 

% Simulation time
param.phase = 0;           % current phase
param.nt = 6;              % number of phases
param.dT = [5 5 5 5 5 5];  % duration of each phase
rng(101)
format long

% Observation and forward operator
param.h = @(sol) extract_obs_tr(sol);
param.f = @callTOUGH2; 

% numerical differentiation parameters
param.dh = 1e-2.*[1e7;1;29];
param.df = 1e-4.*[1e7;1;29];

% functions to initialize error structure
param.R = @exp5_getR;
param.initCOV = @COV_comp_FMM;

% Initialization
param.pm0 = 1e-13; % base value of permeability, see INFILE
param.icfile = [obs_path,'/tracer_ic.mat'];
a = load(param.icfile);
iSol = a.iSol;  

iSol.transform = @transform_sol_tr;
iSol.vec = iSol.transform(iSol,param);
