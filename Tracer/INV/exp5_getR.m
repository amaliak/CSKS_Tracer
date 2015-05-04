function R = exp5_getR(param)
vec = ones(param.nmeas,1);
for i = 1:length(param.obsindex)
    vec(param.obsindex{i}) = param.obsstd(i).^2.*ones(length(param.obsindex{i}),1);
end
R = diag(vec);
end
