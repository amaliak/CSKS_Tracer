%function cn=exp1inv(scalep)
% MAIN FILE FOR COMPRESSED STATE KALMAN SMOOTHER FOR LARGE SCALE PROBLEMS
% AMALIA KOKKINAKI, JUDITH YUE LI, PETER KITANIDIS,  STANFORD UNIVERSITY
% 
% Modified for inversion of tracer test
% Added to git (05/03/2015)

clear all;
%dbstop if error
setCSKFpath;
time_KF = tic;
param.rank = 30
disp([datestr(now),'Initialization started ... ']);
foldername = 'TestingNewIO_Rank30_test'
mkdir(foldername);
disp([datestr(now),'Case: ',foldername]);

%%%%%%% INITIALIZATION %%%%%%%%% 
% set to 1 if running from the beginning and need to regenerate bases
reinit = 1; 

time_init = tic; 
if reinit == 1
	% Setting up observations and true solution
	 
	load('../INIT/truemodelsol1_tr.mat','tSols','tObs','modpmxun2')
	%pmx_initial = modpmxun2;

	nmeas = length(tObs{2});
	obsindex = {1:40;41:80;81:120};

	[param,iSol]=modelsetting(nmeas,obsindex,obs_path,param);

	% Initialize covariances and bases
	[U,Sigma,V,~,~,~] = param.initCOV(param); 

	rng(100);
	Sol = iSol ; 
	Solx = 1.*ones(param.m,1);
	%Solx = pmx_initial; 
	stat.pmerr0 = norm(Solx - tSols{1}.pmx)/norm(tSols{1}.pmx);
	disp('Initial permeability error is:')
	disp(stat.pmerr0)
	Sol.pmx = embed_KF_sol(Solx,Sol.pmx,param.x,param.y,param.z,1);
	Sol.vec = Sol.transform(Sol,param);
	
	% Obtain predicted measurements at k = 1
	
	fobs = TOUGH2obs(Sol,param); 
	
	save(['exp3inv_init_r_',num2str(param.rank),'.mat'])
else
	load(['exp3inv_init_r_',num2str(param.rank),'.mat'])
end

% parameters that are frequently changed between runs and cannot be loaded

% UNSCALED STANDARD DEVIATION OF OBSERVATIONS, R matrix

param.obsstd=[1000;0.2;1]; % for unscaled values Case 5

% SCALING

param.x_scale = param.obsstd;

% SCALED STANDARD DEVIATION OF OBSERVATIONS
param.obsstd=abs(param.obsstd./param.x_scale);

R = param.R(param); % identity

% Regularization
param.theta=[1 1];
Sigma = param.theta(1)*param.theta(2)*Sigma;
R = param.theta(2) * R; 
time_init=toc(time_init);
disp([datestr(now),'Time for initialization ',num2str(time_init),' seconds']);

%%% KALMAN FILTER STARTING %%%%%%%%%%%%%%%%
close all
allsteps = 1; 

for k=1:2 %param.nt
	time_1step = tic; 
	disp([datestr(now),' Step K = ',num2str(k),' started.'])

	param.phase = k; 
	savefilename = ['./',foldername,'/sol_step',num2str(k),'.mat'];

	disp([datestr(now),'... HU COMPUTATION Step ..',num2str(k)])

	time_HU = tic;
	param.recalcHU=1;  
% Compute Jacobian around best forecast for Sol(t=0) using fobs(t=50days)
	if k==1
	   HU = JacobTOUGH2obs1st(Sol,fobs,U,param); 
	else
	   HU = JacobTOUGH2obs(Sol,fobs,U,param); 
	end
	
	time_HU = toc(time_HU);
	disp(['Finished computing HU in ',num2str(time_HU),'seconds']);
	
% Scale the observations, contaminate the observations with error
	
	[y,Hx,HU] = ScaleObs(tObs{k+1},fobs,HU,param);

% Analysis = update state at t=k with info from  measurements at t=k+1
% Gives Sol(t=k-1)), Sigma(t=k-1)
	disp([datestr(now),'...ANALYSIS/SMOOTHING Step .. ',num2str(k)])
% the only place where the logk is changed is in the analysis		
	[Sol,Sigma,stat,KG,cn] = FKFanalysis(Sol,Sigma,HU,U,R,Hx,y,param);
	% FKF analysis updates the whole Sol structure for step k-1, including parameters and states. The analyzed states stored in structure Sol, are then used in TOUGH2update as the INCON
	% At K=2, Sol holds the state
	disp(['Step',num2str(k),' after analysis ',num2str(length(Sol.X(Sol.X<0))),' are below zero']);
	Sol.X(Sol.X<0)=0;
	disp(['Step',num2str(k),' after analysis ',num2str(length(Sol.X(Sol.X>1))),' are above one']);
	Sol.X(Sol.X>1)=1; 
    	movefile('K.mat',['./',foldername,'/KG_',num2str(k),'.mat'])

	Solx = extract(Sol,6,param.x,param.y,param.z,1);
	stat.perr = norm(Solx.pressure-tSols{k}.pressure)/norm(tSols{k}.pressure);
    	stat.cerr = norm(Solx.X-tSols{k}.X)/norm(tSols{k}.X);
	stat.pmerr = norm(Solx.pmx - tSols{k+1}.pmx)/norm(tSols{k+1}.pmx); 
	
% Sol0 is s_k-1|k
	Sol0 = Sol; 

% Now propagate Sol to current time k by running TOUGH2 and using s_k-1|k for i.c.
	disp([datestr(now),' FORECAST Step ..  ',num2str(k)]);
	Sol = TOUGH2update(Sol,param); 
	fSol = Sol; % s_k+1|k+1
	plotX1C(fSol,2,['Solution_Step_',num2str(k)]);
	movefile(['Solution_Step_',num2str(k),'.png'],['./',foldername]); 	
	
% check measurements reproduction and plot residuals
	cobs = param.h(fSol); % param.h only extracts observations does not run TOUGH2
	cobs.vec = [cobs.pressure;cobs.X; log(cobs.pmx)+log(param.pm0)]; 
	trobs = y ;% scaled data with noise
        plotres(cobs,trobs,fobs,foldername,obsindex,param);
	
% Calculate FU for  to calculate posterior covariance	
	disp([datestr(now),'... FU CALCULATION Step ..',num2str(k)])
	
	time_FU = tic;
% FU calculation done for current time t=k, Sol0=Sol(t=k-1) used as initial condition
	recalcFU =1 ; 

	if recalcFU == 1
		FU =JacobTOUGH2(Sol0,Sol,U,param);
		save('FU.mat','FU');
	else
		load('FU.mat','FU');
	end
	time_FU = toc(time_FU);
	disp(['Finished computing FU in ',num2str(time_FU),' seconds']);

	disp([datestr(now),'GETTING POSTERIOR COVARIANCE Step ..',num2str(k)])
% Propagate uncertainty Sigma(t=k) = Sigma(t=k-1) + R ..
	Sigma = FKFforecast(Sigma,FU,U,V);
	fSigma = Sigma;

% norm of posterior covariance should be decreasing gradually
	stat.nsigmap = norm(fSigma(1:param.r,1:param.r));
	stat.nsigmac = norm(fSigma(param.r+1:2*param.r,end));	
    	stat.nsigmak = norm(fSigma(2*param.r+1:3*param.r,end));	

	save(savefilename)
	
	if k < param.nt
		disp('%%%%%%...... PREPARATION FOR NEXT STEP ....%%%%%%%')
		disp([datestr(now),'Forecast observation at step',num2str(k+1)]);
		fobs = TOUGH2obs(Sol,param); 
		time_1step = toc(time_1step);
	end

	%eval(['stat_',num2str(k),'=stat;']);
	statistics{k}=stat;
	
	%if k ==1 
	%	save(['./',foldername,'/statistics.mat'],['stat_',num2str(k)]);
	%else
	%	save(['./',foldername,'/statistics.mat'],['stat_',num2str(k)],'-append');
	%end

	save

	disp([datestr(now),' Step ',num2str(k),' finished.'])
	disp('--------------------------------------------------')
	disp('--------------------------------------------------')
end

save(['./',foldername,'/statistics.mat'],'statistics');
% plot statistics
if allsteps
cd(foldername)
plotstats(statistics)
cd ..
end

time_KF = toc(time_KF)

disp([datestr(now),'All ',num2str(k),' steps finished'])
