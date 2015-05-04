function plotwells()

zmax=1019;%-1019;
zmin=1001;%-1019;

[Xiw,Yiw,Ziw]=cylinder([1 1],10);
Ziw=Ziw+zmin;
Ziw(2,:)=zmax+2;
surf(Xiw+321,Yiw+321,Ziw,'FaceColor','r')

% other injection wells
iw=[351 321;291 321;321 291;321 351];
for iwell=1:length(iw)
    surf(Xiw+iw(iwell,1),Yiw+iw(iwell,2),Ziw,'FaceColor','r')
end

% monitoring wells
mw=[291 291; 291 351; 351 291; 351 351; 381 321; 261 321; 321 381; 321 261];
for iwell=1:length(mw)
    surf(Xiw+mw(iwell,1),Yiw+mw(iwell,2),Ziw,'FaceColor','b')
end