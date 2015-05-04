% function that adds all paths containing CSKF related m-files into the current path, so that all functions can be used.

%% Specify PATH
Main_path = '/home/amaliak/CSKFtracer';
obs_path = '/home/amaliak/Tracer/INV'; % store getObs, getR subroutine
mTOUGH_path = [Main_path,'/mTOUGH'];
kfTOUGH_path = [Main_path,'/kTOUGH'];
mexBBFMM3D_path = [Main_path,'/mexBBFMM3DU'];
addpath(obs_path,mTOUGH_path,kfTOUGH_path,mexBBFMM3D_path);
