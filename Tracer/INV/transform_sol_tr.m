function output = transform_sol_tr(varargin)

% transform a TOUGH2 struct in to a vector for inversion
% Usage: vec = transform(Sol,param);
% Feb 26, 2015 adjusted for tracer

if nargin == 2 && isstruct(varargin{1}) == 1
    
    eSol = varargin{1};
    param = varargin{2};
    p_T = eSol.pressure;
    %temp_T = eSol.T; 
    %s_T = erfinv(2*eSol.s(:,1)-1);
    X2_T = eSol.X;
    pmx_T = log(param.pm0) + log(eSol.pmx);
    output = [p_T;X2_T;pmx_T];

    
% update the TOUGH2 struct to the value of a vector 
% Usage: Sol = transform(vec,Sol,param);
elseif nargin == 3 && isnumeric(varargin{1}) == 1 && isstruct(varargin{2})
    % back transform to x
    input1 = varargin{1};
    eSol = varargin{2};
    param = varargin{3};
    M = param.M;
    eSol.vec = input1;
    eSol.pressure = input1(1:M);
    eSol.X = input1(M+1:2*M);
    %s_T = input1(m+1:2*m);
    %sat = 0.5*(erf(s_T)+1);
    %eSol.s = [sat 1-sat];
    %temp_T = input1(2*m+1:3*m);
    pm_T = input1(2*M+1:end);
    eSol.pmx = exp(pm_T)./param.pm0;
    output = eSol;
else
    disp('Input for transform_sol is in wrong format');
    keyboard;
end

end
