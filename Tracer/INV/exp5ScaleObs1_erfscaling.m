function [y,Hx,HU,R,R0] = exp5ScaleObs1(tobs,fobs,HU,param)
y0 = tobs.vec;
Hx = fobs.vec;
st2s = @(x) 0.5*(erf(x)+1);
% dstds = @(x) sqrt(pi)*exp(erfinv((2*x-1).^2));
diagR0 = ones(param.n,1);
diagR  = diagR0;
%scale obs
y = y0;

figure;

for i = 1:length(param.obsindex)
    % perturb true observation
    if i==3 %saturation is scaled differently
        obs0 = s2st(tobs.sat);
        % non-uniform true STD in saturation measurements
        tsigma = param.tobsstd(i).*dstds(tobs.sat); 
        disp(num2str(mean(tsigma)));
        obs = s2st(tobs.sat) + tsigma.*randn(size(tobs.sat));
        % scaling using nonuniform R
        sigma = param.obsstd(i).*dstds(tobs.sat);  
        y(param.obsindex{i}) = obs./sigma;
        Hx(param.obsindex{i}) = Hx(param.obsindex{i})./sigma;
        HU(param.obsindex{i},:) = bsxfun(@rdivide,HU(param.obsindex{i},:),sigma);
        diagR0(param.obsindex{i}) = sigma.^2;        
    else
        obs0 = y0(param.obsindex{i});
        obs = obs0 + param.tobsstd(i)*randn(size(obs0));
        % scale by sqrt(R)
        diagR0(param.obsindex{i}) = param.obsstd(i).^2.*ones(length(param.obsindex{i}),1);
        y(param.obsindex{i}) = obs./param.obsstd(i);
        Hx(param.obsindex{i}) = Hx(param.obsindex{i})./param.obsstd(i);
        HU(param.obsindex{i},:) = HU(param.obsindex{i},:)./param.obsstd(i);
    end
    % plot
    subplot(length(param.obsindex),1,i)
    if i==3
        plot(param.obsindex{i},st2s(obs0),'r*', param.obsindex{i},st2s(obs),'o')
    else
        plot(param.obsindex{i},obs0,'r*', param.obsindex{i},obs,'o')
    end
    
end
    R0 = diag(diagR0);
    R  = eye(size(R0));

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
    
    function dx = dstds(x)
        % replace x <1e-3 to 1e-3
        x = x.*(x>1e-3) + 1e-3.*(x<1e-3);
        dx = sqrt(pi)*exp(erfinv((2*x-1).^2));
    end

end