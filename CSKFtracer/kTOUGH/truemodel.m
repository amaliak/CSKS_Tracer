function [tSols,tObs] = truemodel(tSol,param)

% true model in folder CSKFtherm
% Run TOUGH2 to simulate true state and observation for one time step
% produces data (without noise)  and true solution
tSols = cell(1,param.nt);
tObs = cell(1,param.nt);
i = 1;
t = 0;
nt = param.nt;    
dT = param.dT; 
T=sum(dT(1:nt)); % total experiment (physical) time 

%% Assign initial condition

% EXTRACT small domain
% tSol is the full solution
% extracting first 5 fields of initial conditions (P,T,elem,por,pmx)
% tSolx=extract(tSol,5,281:2:359,281:2:359,1019:-2:1001);
tSolx=extract(tSol,6,param.x,param.y,param.z,1);

tSolx.transform = @transform_sol_tr;
tSolx.vec = tSolx.transform(tSolx,param);

% tSols stores the solution that is used in the KF (the extracted domain)
tSols{1} = tSolx; 

%CH = getCH(param.x,param.y,param.z,param.L,param.kernel);

param.init = 1;

while t < T,
    disp(['Period ',num2str(i),'  has started']);
    i = i + 1;
    t = t + dT(i-1);
    %% Forward
    % phase 0 corresponds to time 0-->1st data assimilation point
    
    % next line gives full solution, no extraction
    % tSol is the initial tSol matrix and is only used for the pmx vector
    % tSol1 will be the updated solution, full matrix
    tSol = TOUGH2update(tSol,param);
%    tSol0 = tSol;%%
%    % Add model noise
%    tSol.pressure = tSol.pressure + param.dp_std.*CH*randn(size(tSol.pressure));
%    sat = s2x(tSol.s(:,1)) + param.ds_std.*CH*randn(size(tSol.s(:,1))); sat = x2s(sat);
%    tSol.s = [sat 1-sat];
%    tSol1 = tSol; %%
    
    % extract solution for smaller domain to be used for inversion
    tSolx = extract(tSol,6,param.x,param.y,param.z,1);
    tSolx.transform = @transform_sol_tr; %for scaling
    tSolx.vec = tSolx.transform(tSolx,param);
    
    % tSols stores the solution that is used in the KF (the extracted domain)

    tSols{i} = tSolx;
    
    %% Observation
    obs = TOUGH2obs(tSol,param);
    
    % %% Add data noise
    % obs.flux = obs.flux + param.obsflux_std*randn;
    % obs.pressure = obs.pressure + param.obspres_std*randn;
    
    tObs{i} = obs.vec;
    disp(['STEP ',num2str(param.phase),' DONE'])
    param.phase = param.phase + 1 ; 
end

save('truemodelsol1_tr.mat','tSols','param','tObs');

    function CH = getCH(x,y,z,corlengths,kernelfun)
        %x,y,z correspond to (x,y,z) triplets!
        d=cov_irg(x,y,z,corlengths(1),corlengths(2),corlengths(3));
        Q = kernelfun(d);
        CH = chol(Q + 1e-12*eye(size(Q)))';% force Q to be positive definite
    end
end
