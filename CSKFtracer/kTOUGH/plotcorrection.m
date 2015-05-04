function plotcorrection(KG,dres,state_type,xyplane,step)

% function to plot the correction to a particular state 
% for each type of measurement and cumulatively
% example usage 

% plotKG3var(K,40,'logk','logk',1);
% plotKG3var(K,40,'pressure','logk',1);
% plotKG3var(K,40,'concentration','logk',1);


% mi : number of unknown state of each type (24040)
% ni : number of measurement of each type (e.g. 40)

if strcmp(state_type,'pressure')
    i = 1; 
elseif strcmp(state_type,'concentration')
    i = 2;
elseif strcmp(state_type,'logk')
    i = 3;  
end

figuretitle=['Kalman correction for ',state_type,' at step ',num2str(step)];
filename=['correction_',state_type,'_Step',num2str(step)];
%%%% 
figure

load MESHmap_tSols_pump_small.mat xt yt zt
mi = length(xt);

% extract the correct block of the KG matrix that corresponds to the state
% of interest

KGij = KG(mi*(i-1)+1:i*mi,:); % size mxn
obsindex = [40 40 40];
obstype ={'pressure','concentration','logk'}; 
dS = KGij * dres;
[X,Y,~,dStotal]=vec2mat(xt,yt,zt,dS);

subplot(4,1,1)
contourf(X(:,:,xyplane),Y(:,:,xyplane),dStotal(:,:,xyplane),10,'LineStyle','None')
title(['Total correction in ',state_type])
colorbar
cmin = min(min(dStotal(:,:,xyplane)));
cmax = max(max(dStotal(:,:,xyplane)));
daspect([1 1 1])

for j = 1 : length(obsindex)
    ni = obsindex(j);
    dresij = dres(ni*(j-1)+1:j*ni);
    corij = KGij(:,ni*(j-1)+1:j*ni) * dresij;
    [X,Y,~,dS]=vec2mat(xt,yt,zt,corij);
    subplot(4,1,j+1)
    contourf(X(:,:,xyplane),Y(:,:,xyplane),dS(:,:,xyplane),10,'LineStyle','None')
    dStotal = dStotal + dS;


    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')

    hold on

    pw=[351 351; 351 291; 291 351; 291 291];
    kw=[378 378; 378 262; 262 378; 262 262];
    cw=[319 321; 349 321; 321 289; 321 349];

    if strcmp(obstype{j},'pressure')
        plot(pw(:,1),pw(:,2),'s','MarkerFaceColor','m','MarkerEdgeColor','y')
    elseif strcmp(obstype{j},'concentration')
        plot(cw(:,1),cw(:,2),'s','MarkerFaceColor','g','MarkerEdgeColor','y')
    elseif strcmp(obstype{j},'logk')
        plot(kw(:,1),kw(:,2),'s','MarkerFaceColor','c','MarkerEdgeColor','y')
    end
    % injection well
    plot(321,321,'v','MarkerFaceColor','y','MarkerEdgeColor','b')
    title(['Correction in ',state_type,' due to ',obstype{j}])
    colorbar
    caxis([cmin cmax])
    daspect([1 1 1])

end


% create title
casename = 'Case1Rank30';
foldername = ['~/Tracer/INV/',casename,'/'];
savefig([foldername,filename])
