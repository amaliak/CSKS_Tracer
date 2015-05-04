function plotobs(fobs,filename)

subplot(3,1,1)
plot(1:10,fobs(1:10),'ro',11:20,fobs(11:20),'bo',21:30,fobs(21:30),'go',31:40,fobs(31:40),'mo')
legend('pressure1','pressure2','pressure3','pressure4');
ylabel('Pressure (kPA)')

subplot(3,1,2)
plot(1:10,fobs(41:50),'ro',11:20,fobs(51:60),'bo',21:30,fobs(61:70),'go',31:40,fobs(71:80),'mo')
legend('conc1','conc2','conc3','conc4');
ylabel('C/C_0')

subplot(3,1,3)
plot(1:10,fobs(81:90),'ro',11:20,fobs(91:100),'bo',21:30,fobs(101:110),'go',31:40,fobs(111:120),'mo')
legend('pmx1','pmx2','pmx3','pmx4');
ylabel('logk')


print('-dpng',filename)

end
