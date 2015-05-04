function [Sol,Sigma,stats]=FKFupdate(Sol1,Hx,y,U,Sigma,V,R,FU,HU,param)
% OUTPUT: nonscaled matrix-vector product FU and HU using finite difference scheme
% FU: dynamic sensitivity matrix
% HU: observation sensitivity matrix

%U: preselected basis
%tObs: data
%s: state vector
%delta: a vectot of perturbation coefficients
%Sol0:  s(t|t)
%Sol1 : s(t+1|t) best prediction

%% update prior covariance P = ACA' from P = U Sigma U'
A = U'*(FU);    % mN^2
Sigma = A*Sigma*A' + V; % 2N^3
% innovation
% obs = TOUGH2obs(Sol1);
% y = [tObs.pressure; tObs.flux]; Hx = [obs.pressure;obs.flux];
innovation = y - Hx;

% Kalman gain
HPHT = HU*Sigma*HU';
PSI = HPHT + R; % 2nN^2, innovation matrix
US = U*Sigma; %mN^2
PHT = US*HU';
K = US*HU'/PSI; % n^3 + mNn

% measurement update
Sol = updateSol(Sol1,K*innovation);
Sigma = (eye(size(Sigma))-U'*K*HU)*(US'*U); % Nmn + nN^2 + mN^2 + N^3
% sprintf('the delta is %d', delta(2))
% sprintf('the rank of HA is %d',rank(HA))
a = PSI\innovation;
stats.J = 0.5*a'*HPHT*a + 0.5*innovation'/R*innovation; % -log(posterior)
stats.CR = 0.5*log(trace(PSI)) + 0.5*innovation'*a + 0.5*(2*param.m)*log(2*pi);% -log(p(res)) should be small
stats.Q2 = innovation'*a/param.n; % should be close to 1
stats.res = innovation;
stats.mres = mean(innovation);
stats.stdres = std(innovation);
stats.corrmtx = corrcoef([innovation(2:end) innovation(1:end-1)]);
    function Sol = updateSol(Sol,ds)
        s0= [Sol.pressure;erfinv(2*Sol.s(:,1)-1)];
        s0 = s0 + ds;
        if sum(isnan(s0))>0
            error('State variables contain NaN');
        end
        Sol.pressure = s0(1:end/2);
        sat = s0(end/2+1:end);
        sat = 0.5*(erf(sat)+1);
        Sol.s = [sat 1-sat];
    end

end
