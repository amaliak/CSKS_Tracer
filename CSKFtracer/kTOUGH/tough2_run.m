% Specify PATH
Main_path = pwd;
% mexBBFMM2D_path = '/home/judith/Package/mexBBFMM2D';
dct2d_path = 'DCT2D/';
addpath(Main_path, dct2d_path);
%% Set up folder
load('start.mat')
param.r = 20;
% creating_directories(4*param.r); % 4*r
% % removing_directories(80);
% %% Run TOUGH2 true case
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% produce data and true solution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nt = 12;  %12 x 90 days
% dT = 7776000;%90*day();
% T = nt*dT;
% tSols = cell(1,nt);
% tObs = cell(1,nt);
% i = 0;
% t = 0;
% while t < T,
%     %% Forward simulation
%     cd tmp/
%     if t==0
%         delete INCON
%         system('mpirun -np 10 tough2-mp-eco2n');
%         [dumb,tSol] = dumbSAVE();
%     else
%         tSol = TOUGH2update(dumb,tSol);
%     end
%     cd ../    
%     plotX(tSol);
%     obs = TOUGH2obs(dumb,tSol);
%     %     %% Add data noise
%     %     obs.presvar = (2*barsa())^2;
%     %     obs.fluxvar = (1*meter^3/day())^2;
%     i = i + 1;
%     t = t + dT;
%     tSols{i} = tSol;
%     tObs{i} = obs;
% end
% data.obs = tObs;
% data.nt = nt;
% data.dT = dT;
% data.T = T;
% save('start.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define kernel for Q
% syms h; kernel = exp(-sqrt(h.^0.5)./30); %input kernel as a symbolic function
% getQH = makeFMM2D(mexBBFMM2D_path);
%Initialize QV
delta = [1e-4;1e-2];
% r = 30
param.L = 100;
param.n = 13; param.ny = 45; param.nx = 45; param.m = param.nx*param.ny;
flux_std = 3e-7; flux_scale = 6e-4; fluxvar = (flux_std/flux_scale)^2;
pres_std = 7e1; pres_scale = 6e7; presvar = (pres_std/pres_scale)^2;
param.var = [presvar*ones(9,1);fluxvar*ones(4,1)];
param.pres_scale = pres_scale; param.flux_scale = flux_scale;
rng(100)

param.theta = 2.8e-4;
param.var = param.var;
%% P0, and R
[U,Sigma,R] = init_mrst_noQ(param);
% delta = [1e-3;1e-2];

%% Main Loop
Sigma = param.theta.*Sigma;
R = bsxfun(@times,param.var,R);

% start with a wrong saturation 0.3
t = 0; T = data.T; nt = data.nt; dT = data.dT;
sol = cell(nt,1);
m = 2025;
sat = 0.3*ones(m,1); pressure =  20700000*ones(m,1);
Sol.s = [sat 1-sat]; Sol.pressure = pressure;
%%
for i = 1:nt
    % observation
    obs = data.obs{i}; % change to observation
    obs.flux = obs.flux + flux_std*randn;
    obs.pressure = obs.pressure + pres_std*randn;
    % forecast
    Sol0 = Sol;
    Sol1 = TOUGH2update(dumb,Sol0); fobs = TOUGH2obs(dumb,Sol1); plotX(Sol1)
    res_f = [(fobs.pressure-obs.pressure)./param.pres_scale;(fobs.flux-obs.flux)./param.flux_scale];
    J0 = 1e20;
    for j = 1:1
        % correction
        [Sol,Sigma,U]=FKFupdate_noQ(Sol0,Sol1,obs,U,Sigma,R,delta,dumb,param);
        % compute residual
        eSol = norm(Sol.s(:,1)-Sol1.s(:,1))/norm(Sol1.s(:,1));
        obs_p = TOUGH2obs(dumb,Sol);
        sol{i} = Sol;
        %% Error analysis
        load stats.mat
        pobs = TOUGH2obs(dumb,Sol);
        res = [(pobs.pressure-obs.pressure)./param.pres_scale;(pobs.flux-obs.flux)./param.flux_scale];
        Pe = HPHT+R; Peinv = inv(Pe); a = Peinv*res;
        J = 0.5*a'*HPHT*a + 0.5*res'*inv(R)*res; % -log(posterior)
        if J > J0 || eSol < 1e-3
            disp(['converge in ',num2str(j),' runs'])
            if J > J0
                Sol = Sol1;
            end
            break;
        else
            J0 = J;
            Sol1 = Sol;
        end
    end
    
    C = 0.5*log(trace(Pe)) + 0.5*res'*a + 0.5*param.m*log(2*pi);% -log(p(res)) should be small
    Q2 = res'*a/param.n; % should be close to 1
    MSEf = norm(res_f);
    MSEp = norm(res)
    Err = norm(Sol.s(:,1)-tSols{i}.s(:,1))/norm(tSols{i}.s(:,1));
    %% Plot
    plotX(Sol1) % plot the corrected image
    plotX(tSols{i}) % plot the true image

end



