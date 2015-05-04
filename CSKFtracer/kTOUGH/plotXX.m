function plotXX(tSols)

[X,Y,Z,K]=vec2mat(tSols{1}.pmx);

% transform Z

for i=2:6
    
figure;
[X,Y,Z,P]=vec2mat(tSols{i}.pressure-10^7);

h=slice(X,Y,Z,K*10^4,[],[],-18);
set(h,'LineStyle','none');
hold on
h2=slice(X,Y,Z,P,321,321,[]);
set(h2,'LineStyle','none');

plotwells()

hold off
axis tight
daspect([1 1 1])
view(gca,[39.5 26]);
print('-dpng',['threeD_',num2str(i)])

end
