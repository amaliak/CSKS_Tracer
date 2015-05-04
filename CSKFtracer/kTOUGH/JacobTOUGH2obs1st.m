function [HU]=JacobTOUGH2obs1st(Sol1,obs,U,param)
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
%Sol0:  s(t|t)
%Sol1 : s(t+1|t) best prediction

HU = zeros(param.n,sum(param.r));
recalcHU=param.recalcHU;
%% Variable 3
col = [sum(param.r(1:2))+1:sum(param.r)];
if abs(param.dh(3))>1e-6 && recalcHU
    HU3 = observation1(Sol1,obs,U(:,col),param.dh(3)); % or from observation
    save('HU3.mat','HU3');
    HU(:,col) = HU3.forward.vec;
elseif ~recalcHU
    load('HU3.mat','HU3');
    HU(:,col) = HU3.forward.vec;
else
    fprintf('dh(3) is too small\n');
    keyboard;
end
% sprintf('the delta is %d', delta(2))
% sprintf('the rank of HA is %d',rank(HA))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%function%%%%%%%%
    function Sol_perturbed = updateSol(Sol,dsx,param)
	% Input: Sol is the current estimate of the unknown state and parameters FULL domain
	% Output: Sol is the perturbed estimate of the state and parameters FULL domain
	% dsx is the perturbation, x: corresponds to the extracted domain
	s0 = Sol.vec;
	M  = param.M;
	m = param.m;
	% s0 is 2 x (1:M) (M: bid domain)
	% ds is 2 x (1:m) ds(dpressure,dpmx)

    % ds for first and second variable is zeros
    %ds(1:M,1) = embed_KF_sol(dsx(1:param.m),zeros(M,1),param.x,param.y,param.z,1);
    %ds(M+1:2*M,1) = embed_KF_sol(dsx(m+1:2*m),zeros(M,1),param.x,param.y,param.z,1);
    ds(2*M+1:3*M,1) = embed_KF_sol(dsx(2*m+1:3*m),zeros(M,1),param.x,param.y,param.z,1);
    ds(1:M,1) = zeros(size(ds(2*M+1:3*M,1)));
    ds(M+1:2*M,1) = zeros(size(ds(2*M+1:3*M,1)));

    s0 = s0 + ds;
    if sum(isnan(s0))>0
        disp('there is NAN is the state vector')
        keyboard;
    end

    % s0 has the perturbed P and logK that will be used to run TOUGH2
    % Sol is the unperturbed solution as a structure
    % Sol_perturbed is the new perturbed solution as a structure with only P and logK perturbed

    Sol_perturbed = Sol.transform(s0,Sol,param); % takes log(m2) and gives modifier

    end

    function HU = observation1(Sol1,obs,U,a)
        N = size(U,2); 
        n = length(obs.vec);
        Y0= repmat(obs.vec,1,N); % observations we already have for step k based on Sol1, Y0: n x r
        % perturb each Sol1 by delta*a*U
	% perturbation should only happen in the domain to be inverted
	% extraction and embedding will happen within updateSol
        eSol1 = cell(N,1);% eSol2 = cell(N,1);
        for i = 1:N
		disp(['Perturbed realization ',num2str(i)]) 
        	eSol1{i} = updateSol(Sol1,a*U(:,i),param); % perturb Sol1 --> eSol1
%       	eSol2{i} = updateSol(Sol1,-a*U(:,i));
        end
        % simulate N observation in parallel
        % each eSol1 has a perturbed logk 
        %(AK) obs1 corresponds to a propagation of the perturbed eSol1 (initial condition, t=0) in
        %time i.e. the current time &
        %%
	% $x^2+e^{\pi i}$ also obs is at current time
	
        obs1 = TOUGH2obs_parallel(eSol1,param,N); % run TOUGH2 to get obs for eSol1
%       obs2 = TOUGH2obs_parallel(eSol2,param,N);
        Y1 = zeros(size(Y0)); %Y2 = zeros(size(Y0));
        for i = 1:N
        	Y1(:,i) = obs1{i}.vec;
%         	Y2(:,i) = obs2{i}.vec;
        end
        % compute HA
        HU.forward.vec = bsxfun(@rdivide,(Y1-Y0),a);
%         HU.central.vec = bsxfun(@rdivide,(Y1-Y2),2.*a);
    end

end
