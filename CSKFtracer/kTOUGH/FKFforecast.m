function Sigma = FKFforecast(Sigma,FU,U,V)
%% update prior covariance P = ACA' from P = U Sigma U'
A = U'*(FU);    % mN^2
Sigma = A*Sigma*A' + V; % 2N^3
end
