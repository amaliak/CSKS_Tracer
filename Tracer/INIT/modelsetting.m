%% Specify PATH
Main_path = '/home/amaliak';
obs_path = '/home/amaliak/Tracer/INIT'; % store getObs, getR subroutine
mTOUGH_path = [Main_path,'/CSKFtracer/mTOUGH'];
kfTOUGH_path = [Main_path,'/CSKFtracer/kTOUGH'];
%mexBBFMM2D_path = [Main_path,'/mexBBFMM2D'];
%dct2d_path = [Main_path,'/DCT2D'];
addpath(obs_path,mTOUGH_path,kfTOUGH_path);

% COMMENTING OUT WHATEVER IS NOT NEEDED!

%% Setting up folders and MPI for each delta
% Variables estimated: pressure head and logk
%param.r = [1;1;1]; % number of orthorgonal basis for each variable
%nfolder = max(param.r);

nfolder = 1; 
param.nodes = 1; param.ppn = 24;
param.pps = 24; % #processor per simulation
creating_directories(); % remove all folders starting with obs and test
writing_processors_name(param,nfolder); % create state_update and obs_update and add nodefile 'node1'   % needed to call tOUGH2


%% Setup paramters

% Grid information
%param.n = 200; 
%param.ny = 45; param.nx = 45; 
%param.x = [281:2:359]; param.y = [281:2:359]; param.z = [1019:-2:1001];
load coord_htr2.mat
param.x=x_htr2; param.y=y_htr2; param.z=z_htr2;
param.m = length(param.x) ; %29806; %param.nx*param.ny;


% kernel
param.L = [100;100;10];
param.kernel = @(h) exp(-(h.^2));

% Simulation time: 5x10 days
param.phase = 1;
param.nt = 6; %number of timesteps, for forward simulation only need forst timestep
param.dT = [5 5 5 5 5 5]; % in days up to day 12.5
rng(101)
format long

% Observation and forward operator
param.h = @(sol) extract_obs_tr(sol);
param.f = @callTOUGH2; % assume node1 as node file name, need to be fixed


% Initialization --> this is the true solution
% it also has to contain the true heterogeneous field
% ic.mat is the transformed INCON file
param.pm0 = 1e-13; % from INFILE
param.icfile = [obs_path,'/tracer_ic.mat'];
%param.icfile = [obs_path,'/exp6t_ic.mat'];
a = load(param.icfile);
iSol=a.iSol;

iSol.transform = @transform_sol_tr; %for scaling
iSol.vec = iSol.transform(iSol,param);
