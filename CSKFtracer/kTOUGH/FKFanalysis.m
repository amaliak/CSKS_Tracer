function [Sol,Sigma,stats,K,cn] = FKFanalysis(Sol,Sigma,HU,U,R,fobs,tobs,param)
% Delta residual
%dres = tobs.vec - fobs.vec;
dres = tobs - fobs; 
% Kalman gain
HPHT = HU*Sigma*HU';

%PSI = HPHT + R; % 2nN^2, innovation matrix
PSI = HPHT + R ;
disp('The condition number of HPHT + R is')
cn = cond(PSI);
disp(cn)

[UA,SA,VA]=svd(PSI);
V = UA*SA.^0.5;
disp('Aproximation of PSI error')
norm(V*V' - PSI)

%the HPHT is scaled  by the R, so scaling R by the R also gives I 

US = U*Sigma; %mN^2
PHT = US*HU';
K = PHT/PSI; % n^3 + mNn

save K.mat HPHT PSI US PHT K dres

% measurement update
Sol = updateSol(Sol,K*dres);
Sigma = (eye(size(Sigma))-U'*K*HU)*(US'*U); % Nmn + nN^2 + mN^2 + N^3
% sprintf('the delta is %d', delta(2))
% sprintf('the rank of HA is %d',rank(HA))
a = PSI\dres;
a2 = V\dres; %normalized residuals =V'*r
stats.nres=a2;
stats.mres2=mean(a2);
stats.stdres2 = std(a2); 
stats.autocor = autocorr(a2,1);
stats.J = 0.5*a'*HPHT*a + 0.5*dres'/R*dres; % -log(posterior)
stats.CR = 0.5*log(trace(PSI)) + 0.5*dres'*a2 + 0.5*(2*param.m)*log(2*pi);% -log(p(res)) should be small
stats.Q2 = dres'*a2/param.n; % should be close to 1
stats.res = dres;
stats.mres = mean(dres);
stats.stdres = std(dres);
stats.corrmtx = corrcoef([dres(2:end) dres(1:end-1)]);

function Sol = updateSol(Sol,dsx)
        s0= Sol.vec;%[Sol.pressure;erfinv(2*Sol.s(:,1)-1)];i
	% s0 is for the full domain, ds is just for the domain to be inverted
	% expand ds to match size of s0
	M = param.M;
	m = param.m; 
	ds(1:M,1) = embed_KF_sol(dsx(1:m),zeros(M,1),param.x,param.y,param.z,1);
	ds(M+1:2*M,1) = embed_KF_sol(dsx(m+1:2*m),zeros(M,1),param.x,param.y,param.z,1);
	ds(2*M+1:3*M,1) = embed_KF_sol(dsx(2*m+1:3*m),zeros(M,1),param.x,param.y,param.z,1);

        s0 = s0 + ds;
        if sum(isnan(s0))>0
            error('State variables contain NaN');
            keyboard;
        end
        Sol = Sol.transform(s0,Sol,param);
end

end
