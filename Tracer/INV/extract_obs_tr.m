function obs=extract_obs_tr(eSol)

% function to extract the observation matrix tObs from tSol
% analogous to exp6c_obs.m

% tSol has the full domain size

% struc.cur is just a dummy structure being used because function extract
% needs to read in a structure

%%%%%%%%%% MEASUREMENT 1 %%%%%%%%%%%%%%%%
struc1.pressure=eSol.pressure;
z_pt=[1019:-2:1001];

x_pres=[351*ones(size(z_pt)) 351*ones(size(z_pt)) 291*ones(size(z_pt)) 291*ones(size(z_pt))];
y_pres=[351*ones(size(z_pt)) 291*ones(size(z_pt)) 351*ones(size(z_pt)) 291*ones(size(z_pt))];
z_pres=[z_pt z_pt z_pt z_pt];


pres = extract(struc1,1,x_pres,y_pres,z_pres,0);
obs.pressure = pres.pressure;

%%%%%%%%% MEASUREMENT 2 %%%%%%%%%%%%%%%%

struc2.X = eSol.X;
z_temp = [1019:-2:1001];

x_temp = [319*ones(size(z_pt)) 349*ones(size(z_pt)) 321*ones(size(z_pt)) 321*ones(size(z_pt))];
y_temp = [321*ones(size(z_pt)) 321*ones(size(z_pt)) 289*ones(size(z_pt)) 349*ones(size(z_pt))];
z_temp = [z_temp z_temp z_temp z_temp];

temp = extract(struc2,1,x_temp,y_temp,z_temp,0);
obs.X = temp.X;

%%%%%%%%%% MEASUREMENT 3 %%%%%%%%%%%%%%%%%%%%%%

struc3.pmx=eSol.pmx;
z_pmx=1019:-2:1001;

x_pmx=[378*ones(size(z_pmx)) 378*ones(size(z_pmx)) 262*ones(size(z_pmx)) 262*ones(size(z_pmx))];
y_pmx=[378*ones(size(z_pmx)) 262*ones(size(z_pmx)) 378*ones(size(z_pmx)) 262*ones(size(z_pmx))];
z_pmx=[z_pmx z_pmx z_pmx z_pmx];

perm = extract(struc3,1,x_pmx,y_pmx,z_pmx,0);
obs.pmx = perm.pmx;


%%%%%%%% MEASUREMENT VECTOR %%%%%%%%%%%%%
% note that obs.pmx in obs.vec is transformed to logk after param.h is called
obs.vec = [obs.pressure;obs.X;obs.pmx];
