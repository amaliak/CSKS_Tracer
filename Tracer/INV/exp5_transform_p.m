function output = exp5_transform_p(varargin)

% transform a TOUGH2 struct into a vector for inversion
% Usage: vec = transform(Sol,param);
if nargin == 2 && isstruct(varargin{1}) == 1
    % transform to x_T
    eSol = varargin{1};
    param = varargin{2};
    p_T = eSol.pressure;
    %s_T = erfinv(2*eSol.s(:,1)-1);
    pmx_T = log(param.pm0) + log(eSol.pmx); %pmx --> logk
    %output = [p_T;s_T;pm_T];
    output = [p_T;pmx_T];
% takes pmx and gives logk
    
% update the TOUGH2 struct to the value of a vector 
% Usage: Sol = transform(vec,Sol,param);
elseif nargin == 3 && isnumeric(varargin{1}) == 1 && isstruct(varargin{2})
    % back transform to x
    % takes logk and gives pmx
    input1 = varargin{1};
    eSol = varargin{2};
    param = varargin{3};
    M = param.M;
    eSol.vec = input1;
    eSol.pressure = input1(1:M);
    %s_T = input1(m+1:2*m);
    %sat = 0.5*(erf(s_T)+1);
    %eSol.s = [sat 1-sat];
    pm_T = input1(M+1:end); % permeability as logK
    eSol.pmx = exp(pm_T)./param.pm0; % permeability as the TOUGH2 modifier
    output = eSol;
else
    disp('input is in wrong format');
    keyboard;
end

end
