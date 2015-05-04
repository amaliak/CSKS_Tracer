function plotKG(KG,ni,meas_type,state_type,step)

% function to plot the Kalman Gain for a series of measurements
% example usage

% plotKG3var(K,40,'logk','logk',1);
% plotKG3var(K,40,'pressure','logk',1);
% plotKG3var(K,40,'concentration','logk',1);


% mi : number of unknown state of each type (24040)
% ni : number of measurement of each type (e.g. 40)

% Case specific names - change accordingly
if strcmp(meas_type,'pressure')
j = 1;
elseif strcmp(meas_type,'concentration')
j = 2;
elseif strcmp(meas_type,'logk')
j = 3;
end

if strcmp(state_type,'pressure')
i = 1;
elseif strcmp(state_type,'concentration')
i = 2;
elseif strcmp(state_type,'logk')
i = 3;
end

figuretitle=['Kalman Gain for each',meas_type,' measurement for updating ',state_type,' at step',num2str(step)];
filename=['KG_',meas_type,'_',state_type,'_Step',num2str(step)];
%%%%


load MESHmap_tSols_pump_small.mat xt yt zt
mi = length(xt);
figure;
%set(gcf,'units','normalized','outerposition',[0 0 1 1])
%colorbar_max=max(max(KG));
%colorbar_min=min(min(KG));

% extract the correct block of the KG matrix
% e.g. for K12 is the KG for update in variable 1 using measurement 2
% if for n state size and m measurement number we need to extract
% Kij = K(i*(m-1)+1:i*m,j*(n-1)+1:j*n)

KGij = KG(mi*(i-1)+1:i*mi,ni*(j-1)+1:j*ni);

disp(['Min KG:',num2str(min(min(KGij)))])
disp(['Max KG:',num2str(max(max(KGij)))])


for nmeas = 1:ni

subplot(4,10,nmeas)
% convert to regular grid in matrix form for plotting

[X,Y,~,K]=vec2mat(xt,yt,zt,KGij(:,nmeas));

%  plot
xyplane=mod(nmeas,10);
if xyplane == 0 ; xyplane = 10; end;
%xyplane=11-xyplane;
% X(:,:,1) is at 1019; X(:,:,10)is at 1001

contourf(X(:,:,xyplane),Y(:,:,xyplane),K(:,:,xyplane),10,'LineStyle','None')
%caxis([colorbar_min colorbar_max])

if nmeas ==1
colorbar('Location','WestOutside')
end

title(xyplane)
hold on
%plotmeas()
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')

hold on

% emphasize measurement point
if nmeas<=10
meas_id=1;
elseif nmeas<=20
meas_id=2;
elseif nmeas<=30
meas_id=3;
elseif nmeas<=40
meas_id=4;
end

pw=[351 351; 351 291; 291 351; 291 291];
kw=[378 378; 378 262; 262 378; 262 262];
cw=[319 321; 349 321; 321 289; 321 349];

if strcmp(meas_type,'pressure')
plot(pw(meas_id,1),pw(meas_id,2),'s','MarkerFaceColor','m','MarkerEdgeColor','y')
elseif strcmp(meas_type,'concentration')
plot(cw(meas_id,1),cw(meas_id,2),'s','MarkerFaceColor','g','MarkerEdgeColor','y')
elseif strcmp(meas_type,'logk')
plot(kw(meas_id,1),kw(meas_id,2),'s','MarkerFaceColor','c','MarkerEdgeColor','y')
end
% injection well
plot(321,321,'v','MarkerFaceColor','y','MarkerEdgeColor','b')

end
% create title
annotation(gcf,'textbox',...
           [0.286833333333333 0.947269303201506 0.192447916666667 0.0451977401129944],...
           'String',figuretitle,...
           'EdgeColor','none','FontSize',20,'FitBoxToText','on');
%print('-dpng',filename)
casename = 'Case1Rank30';
foldername = ['~/Tracer/INV/',casename,'/'];
savefig([foldername,filename])
