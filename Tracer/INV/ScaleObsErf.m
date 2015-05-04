function [y,Hx,HU] = ScaleObs(tobs,fobs,HU,param)
y = tobs; 
Hx = fobs.vec;

np = param.obsindex{1};
nc = param.obsindex{2};
nk = param.obsindex{3};

% Measyrement type 1
p_scale = param.x_scale(1); 
HU(np,:) = HU(np,:)./(p_scale); 
y(np) = y(np)./p_scale; %true obs
Hx(np,1) = Hx(np,1)./p_scale; %predicted obs

% Measurement type 2
c_scale = param.x_scale(2);
HU(nc,:) = HU(nc,:)./c_scale;
y(nc) = y(nc)./c_scale;
Hx(nc,1) = Hx(nc,1)./c_scale;

% Measurement type 3
logk_scale = param.x_scale(3);
HU(nk,:) = HU(nk,:)./logk_scale;
y(nk)=y(nk)./logk_scale;
Hx(nk,1) = Hx(nk,1)./logk_scale;

%perturbed scaled observations
y0 = y;
y(np) = y(np) + param.obsstd(1)*randn(size(y(np)));
ync = y(nc) + param.obsstd(2)*randn(size(y(nc)));
ync(ync<0)=0;
ync(ync>1)=1;
y(nc)=ync;

rng(99);
y(nk) = y(nk) + param.obsstd(3)*randn(size(y(nk)));
% need to be adding exactly the same error every time

% plot true solution versus current estimate of solution
figure;
% modified 3/23, plot normal values not normalized values
subplot(1,3,1)
plot(np,y(np)*p_scale,'r*',np,y0(np)*p_scale,'o')
set(gca,'FontSize',14)
title('Pressure')
legend('with noise','without noise')

subplot(1,3,2)
plot(nc,y(nc)*c_scale,'r*',nc,y0(nc)*c_scale,'o')
set(gca,'FontSize',14)
title('Concentrations')
%legend('with noise','without noise')

subplot(1,3,3)
plot(nk,y(nk)*logk_scale,'r*',nk,y0(nk)*logk_scale,'o')
set(gca,'FontSize',14)
title('Log permeability')
%legend('with noise','without noise')

print('-dpng',['plotdata',num2str(param.phase)])

end
