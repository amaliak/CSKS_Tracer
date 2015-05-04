function plotX1T(tSols,xyplane,filename)
% tSols is a single structure
% that has pressure, elem, and pmx

% and can correspond to full domain
% the smaller domain will be extracted
% % need mapping of tSols to (x,y,z), available in MESHmap.mat in G.x, G.y,
% G.z, G.elem

if length(tSols.T)>25000
	pt=strtrim(tSols.elem(1:end-1));
	%disp('load MESH');
	if exist('MESHmap_tSols_pump_full.mat') == 2
        	load MESHmap_tSols_pump_full.mat xt yt zt
	else
    		disp('load MESHmap.mat')
		keyboard;
	    for i=1:length(pt)
        	xt(i)=G.x(strcmp(G.elem,pt(i)));
	        yt(i)=G.y(strcmp(G.elem,pt(i)));
	        zt(i)=G.z(strcmp(G.elem,pt(i)));
	    end
   	 save MESHmap_tSols_pump_full.mat xt yt zt
	end
else	
	pt=strtrim(tSols.elem(1:end));
	%disp('load MESH');
	if exist('MESHmap_tSols_pump_small.mat') == 2
        	load MESHmap_tSols_pump_small.mat xt yt zt
	else
    		disp('load MESHmap.mat')
		keyboard;
	    for i=1:length(pt)
        	xt(i)=G.x(strcmp(G.elem,pt(i)));
	        yt(i)=G.y(strcmp(G.elem,pt(i)));
	        zt(i)=G.z(strcmp(G.elem,pt(i)));
	    end
   	 save MESHmap_tSols_pump_small.mat xt yt zt
	end
end
% then convert to regular grid in matrix form for plotting
%[X,Y,~,K]=vec2mat(tSols{1}.pmx);
[X,Y,~,K]=vec2mat(xt,yt,zt,tSols.pmx(1:end));
% and then plot 

figure; 
subplot(2,1,1)
contourf(X(:,:,xyplane),Y(:,:,xyplane),K(:,:,xyplane),10,'LineStyle','None')
daspect([1 1 1])
colorbar
title(['Permeability modifier at Level ',num2str(xyplane)])
hold on
plotmeas()

subplot(2,1,2)
[X,Y,~,T]=vec2mat(xt,yt,zt,tSols.T(1:end));
contourf(X(:,:,xyplane),Y(:,:,xyplane),T(:,:,xyplane),20,'LineStyle','None')
colorbar('FontSize',12)
%caxis([9.4*10^8 9.9*10^8]);
title(['\Temperature'])
daspect([1 1 1])
hold on
plotmeas()


print('-dpng',filename)

% plot measurement points
function plotmeas()

x_pt=[291 291 351 351];
y_pt=[291 351 291 351];


x_kt=[335 335 307 307];
y_kt=[335 307 335 307];

for i=1:length(x_pt)
	plot(x_pt(i),y_pt(i),'kx')
	plot(x_kt(i),y_kt(i),'ko')
end

%legend(' ','Pressures','Permeabilities','Location','southoutside')
