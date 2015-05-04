function plotXXT(tSols)

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



[X,Y,Z,K]=vec2mat(xt,yt,zt,tSols.pmx(1:end-1));

% transform Z
figure;
[X,Y,Z,T]=vec2mat(xt,yt,zt,tSols.T(1:end-1));

h=slice(X,Y,Z,-log10(K*10^(-12)),258,381,1002);
set(h,'LineStyle','none');
colormap('hot')
hold on
freezeColors
h2=slice(X,Y,Z,T,321,321,[]);
set(h2,'LineStyle','none');
colormap('jet')
colorbar

plotwells()

hold off
axis tight
daspect([1 1 1])
view(gca,[39.5 26]);
%print('-dpng',['threeD_T'])
savefig('threeD_T')

end
