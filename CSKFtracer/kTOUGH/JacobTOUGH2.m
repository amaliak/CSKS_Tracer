function FU=JacobTOUGH2(Sol0,Sol1,U,param)
% OUTPUT: nonscaled matrix-vector product FU and HU using finite difference scheme
% For single variable
% FU: dynamic sensitivity matrix, FU = dfdp*Up =
% (f(p+delta*norm(Up))-f(p))/delta
% HU: observation sensitivity matrix
% c: central difference
% f: forward difference

%U: 2m x r (e.g., dfdp x Up)
%tObs: data
%s: state vector
%delta: a vectot of perturbation coefficients
%Sol0:  best prediction
%Sol1 : forecast

%% Variable 1
col = 1:param.r(1);
if abs(param.df(1))>1e-6
    FU1 = forward(Sol1,Sol0,U(:,col),param.df(1)); % or from observation
    save('FU1.mat','FU1');
    FU(:,col) = FU1.forwardx.vec;
else
    fprintf('df(1) is too small\n');
    keyboard;
end

%% Variable 2
col = (param.r(1)+1):(param.r(1)+param.r(2));
if abs(param.dh(2))>1e-6
    FU2 = forward(Sol1,Sol0,U(:,col),param.df(2)); % or from observation
    save('FU2.mat','FU2');
    FU(:,col) = FU2.forwardx.vec;
else
    fprintf('df(2) is too small\n');
    keyboard;
end

%% Variable 3
col = (sum(param.r(1:2))+1):sum(param.r);
if abs(param.dh(3))>1e-6
FU3 = forward(Sol1,Sol0,U(:,col),param.df(3)); % or from observation
save('FU3.mat','FU3');
FU(:,col) = FU3.forwardx.vec;
else
fprintf('df(3) is too small\n');
keyboard;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%function%%%%%%%%
    function FU = forward(Sol1,Sol0,U,a)
        N = size(U,2);
        % f is augmented state vector with pressure and transformed saturation
        Sol1.vec = Sol1.transform(Sol1,param);
        F0 = repmat(Sol1.vec,1,N);% p(t+1|t)
        % perturb each eSols
        eSol1 = cell(N,1);% eSol2 = cell(N,1);
        for i = 1:N 
              disp(['Perturbed realization ',num2str(i)])
	      eSol1{i} = updateSol(Sol0,a*U(:,i),param);
%             eSol2{i} = updateSol(Sol0,-a*U(:,i));
        end
        
        % propagate N eSols to t+dT in parallel
        tic;
        eSol1 = TOUGH2update_parallel(eSol1,N,param);
        t1 = toc;
	disp(['Extra time for incorporateing forward sims to eSols']);
        disp(['TOUGH2update_parallel finish in ',num2str(t1)])
%         tic;
%         eSol2 = TOUGH2update_parallel(eSol2,N);
%         t2 = toc;
%         disp(['TOUGH2update_parallel finish in ',num2str(t2)])
        F1 = zeros(size(F0));
        for i = 1:N
            F1(:,i) = eSol1{i}.vec;
%             F2(:,i) = eSol2{i}.vec;
        end
        % compute FUi
%         FU.central.vec = bsxfun(@rdivide,F1-F2,2.*a);
	M = param.M; 
        FU.forward.pressure = bsxfun(@rdivide,F1(1:M,:)-F0(1:M,:),a);
	FU.forward.pmx = bsxfun(@rdivide,F1(2*M+1:3*M,:)-F0(2*M+1:3*M,:),a);
	FU.forward.X = bsxfun(@rdivide,F1(M+1:2*M,:)-F0(M+1:2*M,:),a);
	% extract smaller domain for FU
	FU.forwardx = extract(FU.forward,3,param.x,param.y,param.z,1);
	FU.forwardx.vec = [FU.forwardx.pressure;FU.forwardx.X; FU.forwardx.pmx];
    end

    function Sol = updateSol(Sol,dsx,param)
        s0= Sol.vec;
	% fixing domain sizes
	M = param.M;
	m = param.m; 
	
	ds(1:M,1) = embed_KF_sol(dsx(1:param.m),zeros(M,1),param.x,param.y,param.z,1);
	ds(M+1:2*M,1) = embed_KF_sol(dsx(m+1:2*m),zeros(M,1),param.x,param.y,param.z,1);
    ds(2*M+1:3*M,1) = embed_KF_sol(dsx(2*m+1:3*m),zeros(M,1),param.x,param.y,param.z,1);

	s0 = s0 + ds;
        if sum(isnan(s0))>0
            disp('there is NAN is the state vector')
            keyboard;
        end
        Sol = Sol.transform(s0,Sol,param);
    end
end
