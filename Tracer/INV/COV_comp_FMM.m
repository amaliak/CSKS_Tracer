function [U,Sigma,V,err,K,CH] = COV_comp_FMM(param)
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

fprintf('# of Gaussian basis for each state variable is %d, %d, %d\n',param.r(1),param.r(2),param.r(3));
% Observation error covariance R
% this is done outside of this function
%R = param.R(param);

%% Initial state error covariance P0
%[x,y] = meshgrid(param.x,param.y);
% kernel = @(h) exp(-h.^2./param.L^2);% correlation length is 7L/4, around 200m
%K = getQ(x(:),y(:),param.kernel);
%CH = chol(K + 1e-14*eye(size(K)))';% force Q to be positive definite

[K,CH]=getCH(param.x,param.y,param.z,param.L,param.kernel);

%% Run FMM to create U,S,V (V not needed)
% FMM needs to be run anew only if the kernel is changed
% and if the number of basis is changed
new_kernel= 1;

if new_kernel ==1
	delete Omega.mat
	[U1,S11,~] = RandomizedCondSVDFMM(param.m,param.r(1));
	[U2,S22,~] = RandomizedCondSVDFMM(param.m,param.r(2));
	[U3,S33,~] = RandomizedCondSVDFMM(param.m,param.r(3));
	delete Omega.mat
	save FMM.mat U1 U2 U3 S11 S22 S33
else
	load FMM.mat U1 U2 S11 S22
end

%% 1. assembling U
%U1 = gen_dct2basis(param.ny,param.nx,param.r(1)); % size mxr
%U2 = gen_dct2basis(param.ny,param.nx,param.r(2));
O2 = zeros(size(U2)); 
O1 = zeros(size(U1));
O3 = zeros(size(U3));
U = [U1 O2 O3;...
     O1 U2 O3;...
     O1 O2 U3]; %2m x (r1+r2)

%% 2. assembling Sigma
S11 = param.x_std(1)^2.*zeros(param.r(1));%no error in initial pressure
%S22 = param.x_std(2)^2.*zeros(param.r(2));%...
%S33 = param.x_std(3)^2.*U3'*K*U3;% use mexBBFMM
S22 = param.x_std(2)^2.*zeros(param.r(2));% no error in initial temperature
S33 = param.x_std(3)^2.*S33;

S12 = param.x_std(1)*param.x_std(2).*zeros(param.r(1),param.r(2));
S13 = param.x_std(1)*param.x_std(3).*zeros(param.r(1),param.r(3));
S23 = param.x_std(2)*param.x_std(3).*zeros(param.r(2),param.r(3));


Sigma = [S11 S12  S13;...
        S12' S22  S23;...
        S13' S23' S33];

%% compressed dynamic error covariance Q
% to make the multiplication U'KU should use FMM again
% need to calculate K*U1 and K*U2 using FMM, so need to supply U1 and U2 to FMM
KU1 = callfmm(U1,0);
KU2 = callfmm(U2,0);
KU3 = callfmm(U3,0);

V11 = param.dx_std(1)^2.*U1'*KU1;%dynamic error in var1 pressure
V22 = param.dx_std(2)^2.*U2'*KU2;%dynamic error in var2 logkdb
V33 = param.dx_std(3)^2.*U3'*KU3;%dynamic error in var3

V12 = param.dx_std(1)*param.dx_std(2).*zeros(param.r(1),param.r(2));
V13 = param.dx_std(1)*param.dx_std(3).*zeros(param.r(1),param.r(3));
V23 = param.dx_std(2)*param.dx_std(3).*zeros(param.r(2),param.r(3));

V = [V11 V12  V13;...
    V12' V22  V23;...
    V13' V23' V33];...

% V = zeros(size(Sigma)); % no dynamic error

% compute relative approximation error;
%O = zeros(param.m,param.m);
%P0 = [ O O;...
%       O param.x_std(2).^2.*K]; % real initial P
P0 = param.x_std(3)^2.*K ;
%P0p = U*Sigma*U'; % approximated initial P
P0p = U3 * S33 * U3' ;
err = norm(P0p-P0,'fro')/norm(P0,'fro');
disp([sprintf('the relative approximation error of P0 is %d',err)]);
varerr = trace(P0p)/trace(P0);
disp([sprintf('%d of the total variance is explained',varerr)]);

function [Q,CH] = getCH(x,y,z,corlengths,kernelfun)
	%x,y,z correspond to (x,y,z) triplets!
        d=cov_irg(x,y,z,corlengths(1),corlengths(2),corlengths(3));
        Q = kernelfun(d);
        CH = chol(Q + 1e-12*eye(size(Q)))'; % forces Q to be positive definite
end

end
