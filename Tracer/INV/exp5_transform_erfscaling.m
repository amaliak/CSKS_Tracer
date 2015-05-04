function output = exp5_transform(varargin)

% transform a TOUGH2 struct in to a vector for inversion
% Usage: vec = transform(Sol,param);
if nargin == 2 && isstruct(varargin{1}) == 1
    % transform to x_T
    eSol = varargin{1};
    param = varargin{2};
    p_T = eSol.pressure;
    % remove any value that >=1 or <=0, and replace it with 1-eps, or eps
    sat = eSol.s(:,1);
    s_T = s2st(sat);
    %     pm_T = log(param.pm0) + log(eSol.pmx);
    output = s_T;%single
    
    % update the TOUGH2 struct to the value of a vector
    % Usage: Sol = transform(vec,Sol,param);
elseif nargin == 3 && isnumeric(varargin{1}) == 1 && isstruct(varargin{2})
    % back transform to x
    input1 = varargin{1};
    eSol = varargin{2};
    param = varargin{3};
    m = param.m;
    eSol.vec = input1;
    s_T = input1;
    sat = st2s(s_T);
    eSol.s = [sat 1-sat];
    %     pm_T = input1(2*m+1:end);
    %     eSol.pmx = exp(pm_T)./param.pm0;
    output = eSol;
else
    disp('input is in wrong format');
    keyboard;
end
    function s_t = s2st(sat)
        for j = 1:length(sat)
            if (sat(j)<=0)
                sat(j) = eps;
            elseif (sat(j) >= 1)
                sat(j) = 1-eps;
            end
        end
        s_t = erfinv(2*sat-1);
    end

    function s = st2s(s_t)
        s = 0.5*(erf(s_t)+1);
    end

end
