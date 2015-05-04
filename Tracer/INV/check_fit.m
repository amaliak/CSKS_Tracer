function cobs=check_fit(Sol,param,y,fobs,name)

% function to check the fit of measurements after smoothing

% Sol:  is the updated solution for the i.c.
% y:    are the contaminated measurements for k+1
% fobs: is the prior prediction s_k+1^f
% name: is used for the image created

Sol = TOUGH2update(Sol,param);

plotX1(Sol,2,['Solution'])

cobs = param.h(Sol);

% check measurements reproduction
cobs.vec = [cobs.pressure; log(cobs.pmx)+log(param.pm0)]; 
trobs = y ;% data with noise
subplot(1,2,1)
plot(1:40,cobs.vec(1:40)/10^7,'r*',1:40,trobs(1:40),'ko',1:40,fobs.vec(1:40)/10^7,'bx')
title('pressures')
legend('predicted','true','prior')
subplot(1,2,2)
plot(1:40,cobs.vec(41:80),'r*',1:40,trobs(41:80),'o')
title('log permeabilities')
legend('predicted','true')
print('-dpng',['residuals',name])

