function plotobs(fobs,filename)

subplot(2,1,1)
plot(1:10,fobs(1:10),'ro',11:20,fobs(11:20),'bo',21:30,fobs(21:30),'go',31:40,fobs(31:40),'mo')
legend('pmx1','pmx2','pmx3','pmx4');
subplot(2,1,2)

plot(1:10,fobs(41:50),'ro',11:20,fobs(51:60),'bo',21:30,fobs(61:70),'go',31:40,fobs(71:80),'mo')
legend('pressure1','pressure2','pressure3','pressure4');

print('-dpng',filename)

end
