function plotres(cobs,trobs,fobs,foldername,obsindex,param)

% this function plots residuals
% cobs are observations after the analysis
% fobs are the observations before the analysis (prior)
% trobs are the true observations

for i=1:length(obsindex)
    eval(['obsindex',num2str(i),'=obsindex{',num2str(i),'};'])
end

figure;
subplot(2,3,1)

plot(obsindex1,cobs.vec(obsindex1),'r*',obsindex1,trobs(obsindex1)*param.x_scale(1),'bo')
hold on
plot(obsindex1,fobs.vec(obsindex1),'r--')
title('Pressures')
legend('updated','true','prior','Position','best')

subplot(2,3,2)
plot(obsindex2,cobs.vec(obsindex2),'r*',obsindex2,trobs(obsindex2)*param.x_scale(2),'bo')
hold on
plot(obsindex2,fobs.vec(obsindex2),'r--')
title('Concentrations')
%legend('predicted','true')

subplot(2,3,3)
plot(obsindex3,cobs.vec(obsindex3),'r*',obsindex3,trobs(obsindex3)*param.x_scale(3),'bo')
hold on
plot(obsindex3,fobs.vec(obsindex3),'r--')
title('log k')
%legend('predicted','true')

%%%%
subplot(2,3,4)
plot(obsindex1,cobs.vec(obsindex1)-trobs(obsindex1)*param.x_scale(1),'x')
title('residuals')
subplot(2,3,5)
plot(obsindex2,cobs.vec(obsindex2)-trobs(obsindex2)*param.x_scale(2),'+')
title('residuals')
subplot(2,3,6)
plot(obsindex3,cobs.vec(obsindex3)-trobs(obsindex3)*param.x_scale(3),'*')
title('residuals')

% save 
print('-dpng','residuals')
movefile(['residuals.png'],['./',foldername,'/residuals_',num2str(param.phase),'.png']);		
movefile(['plotdata',num2str(param.phase),'.png'],['./',foldername,'/plotdata_',num2str(param.phase),'.png']);
