function plotX(tSols,filename)

% % need mapping of tSols to (x,y,z), available in MESHmap.mat in G.x, G.y,
% G.z, G.elem
pt=strtrim(tSols{1}.elem);

if exist('MESHmap_tSols.mat') == 2
    load MESHmap_tSols.mat
else
    load MESHmap.mat
    for i=1:length(pt)
        xt(i)=G.x(strcmp(G.elem,pt(i)));
        yt(i)=G.y(strcmp(G.elem,pt(i)));
        zt(i)=G.z(strcmp(G.elem,pt(i)));
    end
    save MESHmap_tSols.mat xt yt zt
end

% then convert to regular grid in matrix form for plotting
%[X,Y,~,K]=vec2mat(tSols{1}.pmx);
[X,Y,~,K]=vec2mat(xt,yt,zt,tSols{1}.pmx);
% and then plot 

figure; 
subplot(2,3,1)
contourf(X(:,:,5),Y(:,:,5),K(:,:,5),10,'LineStyle','None')
    daspect([1 1 1])

colorbar
title('Permeability modifier')

c=2;
for i=2:6
    subplot(2,3,c)
    [X,Y,~,P]=vec2mat(xt,yt,zt,tSols{i}.pressure-10^7);
    contourf(X(:,:,5),Y(:,:,5),P(:,:,5),20,'LineStyle','None')
    colorbar
    title(['\Delta P at end of phase',num2str(i)])
    daspect([1 1 1])
    c=c+1;
end

print('-dpng',filename)
