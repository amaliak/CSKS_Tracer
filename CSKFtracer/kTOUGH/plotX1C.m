function plotX1(tSols,xyplane,filename)
% tSols is a single structure
% that has pressure, elem, and pmx

% and can correspond to full domain
% the smaller domain will be extracted
% % need mapping of tSols to (x,y,z), available in MESHmap.mat in G.x, G.y,
% G.z, G.elem

% adapted for tracer test inversion

if length(tSols.pressure)>25000
	pt=strtrim(tSols.elem(1:end));
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


[X,Y,~,K]=vec2mat(xt,yt,zt,tSols.pmx(1:end-1));

% and then plot 

figure;
subplot(3,1,1)
contourf(X(:,:,xyplane),Y(:,:,xyplane),log(10^-13*K(:,:,xyplane)),10,'LineStyle','None')
daspect([1 1 1])
colorbar
title('Log Permeability')
hold on
plotmeas()

subplot(3,1,2)
[X,Y,~,P]=vec2mat(xt,yt,zt,tSols.pressure(1:end-1));
%[X,Y,~,P]=vec2mat(xt,yt,zt,tSols.pressure(1:end-1));
contourf(X(:,:,xyplane),Y(:,:,xyplane),P(:,:,xyplane),20,'LineStyle','None')
colorbar
%caxis([9e+06 9.9e+06]);
title(['Pressure (bar) Step'])
daspect([1 1 1])
hold on
plotmeas()

subplot(3,1,3)
[X,Y,~,C]=vec2mat(xt,yt,zt,tSols.X(1:end-1));

%[X,Y,~,P]=vec2mat(xt,yt,zt,tSols.pressure(1:end-1));
contourf(X(:,:,xyplane),Y(:,:,xyplane),C(:,:,xyplane),20,'LineStyle','None')
colorbar
%caxis([9e+06 9.9e+06]);
title(['Concentration [-]'])
daspect([1 1 1])
hold on


plotmeas()




print('-dpng',filename)

% plot measurement points
function plotmeas()

x_pt=[291 291 351 351];
y_pt=[291 351 291 351];

x_tt=[351 321 291 321];
y_tt=[321 351 321 291];

x_kt=[335 335 307 307];
y_kt=[335 307 335 307];

for i=1:length(x_pt)
	plot(x_pt(i),y_pt(i),'kx')
    plot(x_tt(i),y_tt(i),'ks')
	plot(x_kt(i),y_kt(i),'ko')
end

%legend(' ','Pressures','Permeabilities','Location','southoutside')

