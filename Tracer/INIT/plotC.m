function plotC(tSols,level,variable,filename)

% variable can be pressure, T (for temperature) or X (for concentration)



% % need mapping of tSols to (x,y,z), available in MESHmap.mat in G.x, G.y,
% G.z, G.elem
pt=strtrim(tSols.elem(1:end-1));

if exist('MESHmap_tSols_pump_small.mat') == 2
    load MESHmap_tSols_pump_small.mat
else
    disp('MESHmap does not exist')
end

    
figure;
eval(['[X,Y,Z,C]=vec2mat(xt,yt,zt,tSols.',variable,'(1:end));'])

%h=slice(X,Y,Z,-log10(K*10^-12),381,381,1001);
%set(h,'LineStyle','none');
%colormap('hot')
%hold on
%freezeColors
h2=slice(X,Y,Z,C,[],[],level);
%h2=slice(X,Y,Z,-log10(C*10^-12),321,321,1001);

cm=flipud(jet);
cm='jet';
colormap(cm);
%colormap('jet')
colorbar('Location','SouthOutside')
title('Tracer test - Day 10 end of injection')
xlabel(['Level',num2str(level)])
%caxis([10^5 6*10^5])
set(h2,'LineStyle','none');
hold on
plotwells()

hold off
axis tight
daspect([1 1 1])
view(gca,[-54 24]);
print('-dpng',filename)

end
