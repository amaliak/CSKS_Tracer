function plotstats(stats)

nt=length(stats);

for i=1:nt;
    eval(['mres2(i)=stats{',num2str(i),'}.mres2;'])
    eval(['stdres2(i)=stats{',num2str(i),'}.stdres2;'])
    eval(['J(i)=stats{',num2str(i),'}.J;'])
    eval(['Q2(i)=stats{',num2str(i),'}.Q2;'])
    eval(['perr(i)=stats{',num2str(i),'}.perr;'])
    eval(['cerr(i)=stats{',num2str(i),'}.cerr;'])
    eval(['pmerr(i)=stats{',num2str(i),'}.pmerr;'])
    eval(['nsigmap(i)=stats{',num2str(i),'}.nsigmap;'])
    eval(['nsigmak(i)=stats{',num2str(i),'}.nsigmak;'])
    eval(['nsigmac(i)=stats{',num2str(i),'}.nsigmac;'])
end
figure
subplot(2,5,1)
plot(1:nt,mres2,'k*',1:nt,stdres2,'b*','LineStyle','-')
legend('mean','st.dev','Location','NorthOutside')
title('normalized residuals')

subplot(2,5,2)
plot(1:nt,perr,'k*','LineStyle','-')
title('Pressure errors')

subplot(2,5,3)
plot(1:nt,cerr,'k*','LineStyle','-')
title('Conc errors')

subplot(2,5,4)
plot(1:nt,pmerr,'k*','LineStyle','-')
title('Logk errors')

subplot(2,5,5)
plot(1:nt,nsigmap,'k*')
title('Post cov norm for P')

subplot(2,5,6)
plot(1:nt,nsigmac,'b*')
title('Post cov norm for C')

subplot(2,5,7)
plot(1:nt,nsigmak,'r*')
title('Post cov norm for K')

subplot(2,5,8)
plot(1:nt,J,'k*','LineStyle','-')
title('Objective function')

subplot(2,5,9)
plot(1:nt,Q2,'k*','LineStyle','-')
Q2threshold = 1 + 2.8/sqrt(24040-1);
hold on
plot(1:nt,Q2threshold*ones(nt,1),'r--')
title('Q2')

%print('-dpng','Statistics')
savefig('Statistics')

end



