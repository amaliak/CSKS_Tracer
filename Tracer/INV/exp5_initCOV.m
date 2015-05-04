function [U,Sigma,V,err,K,R,CH] = exp5_initCOV(param)
% Initialize error statistics for data assimilation
% P = USU', Q = UVU'
% Qij = theta*kernel(i,j)
% P = zeros(m,m)
% R = var*Inxn
% G: embed geometry
% r: rank
% var: observation variance of  nx1
% theta: regularization parameter of size mx1
%P = [ U2SsU2' U2pSspU2p';
%%%%%  U2pSspU2p' U1SpU1']
%%% where Sigma is full rank 
% var1: s_T
% var2: pressure
% var3: logK

fprintf('# of DCT basis for each state variabe is %d, %d, %d\n',param.r(1),param.r(2),param.r(3));
% Observation error covariance R
R = param.R(param);

%% Initial state error covariance P0
[x,y] = meshgrid(param.x,param.y);
% kernel = @(h) exp(-h.^2./param.L^2);% correlation length is 7L/4, around 200m
K = getQ(x(:),y(:),param.kernel);
CH = chol(K + 1e-14*eye(size(K)))';% force Q to be positive definite

%% 1. assembling U
U1 = gen_dct2basis(param.ny,param.nx,param.r(1)); % size mxr
U2 = gen_dct2basis(param.ny,param.nx,param.r(2));
U3 = gen_dct2basis(param.ny,param.nx,param.r(3));
O2 = zeros(size(U2)); 
O1 = zeros(size(U1));
O3 = zeros(size(U3));
U = [U1 O2 O3;...
     O1 U2 O3;...
     O1 O2 U3]; %3m x (r1+r2+r3)

%% 2. assembling Sigma
S11 = param.x_std(1)^2.*zeros(param.r(1));%no error in initial pressure
S22 = param.x_std(2)^2.*zeros(param.r(2));%...
S33 = param.x_std(3)^2.*U3'*K*U3;% use mexBBFMM
S12 = param.x_std(1)*param.x_std(2).*zeros(param.r(1),param.r(2));
S13 = param.x_std(1)*param.x_std(3).*zeros(param.r(1),param.r(3));
S23 = param.x_std(2)*param.x_std(3).*zeros(param.r(2),param.r(3));
Sigma = [S11 S12  S13;...
        S12' S22  S23;...
        S13' S23' S33];

%% compressed dynamic error covariance Q
V11 = param.dx_std(1)^2.*U1'*K*U1;%dynamic error in var1
V22 = param.dx_std(2)^2.*U2'*K*U2;%dynamic error in var2
V33 = param.dx_std(3)^2.*U3'*K*U3;%dynamic error in var3
V12 = param.dx_std(1)*param.dx_std(2).*zeros(param.r(1),param.r(2));
V13 = param.dx_std(1)*param.dx_std(3).*zeros(param.r(1),param.r(3));
V23 = param.dx_std(2)*param.dx_std(3).*zeros(param.r(2),param.r(3));
V = [V11 V12 V13;...
    V12' V22 V23;...
    V13' V23' V33];
% V = zeros(size(Sigma)); % no dynamic error

% compute relative approximation error;
O = zeros(param.m,param.m);
P0 = [ O O O;...
       O O O;...
       O O param.x_std(3).^2.*K]; % real initial P
P0p = U*Sigma*U'; % approximated initial P
err = norm(P0p-P0,'fro')/norm(P0,'fro');
disp([sprintf('the relative approximation error of P0 is %d',err)]);
varerr = trace(P0p)/trace(P0);
disp([sprintf('%d of the total variance is explained',varerr)]);

	function Q0 = getQ(x,y,kernelfun)
        [xj,xl]=meshgrid(x(:),x(:));
        [yj,yl]=meshgrid(y(:),y(:));
        d = sqrt(((xj-xl)).^2+((yj-yl)).^2);
        Q0 = kernelfun(d);
	end
end