function [Xr,Yr,Zr,Cr]=vec2mat(xt,yt,zt,ct)

% this function takes the 3D irregular grid vector produced by TOUGH2 and 
% creates a 3D matrix of a regular grid 2x2x2 for plotting

%% old code
% x = [281:2:359];
% y = [281:2:359];
% z = [1019:-2:1001];
% z = z - 1019; 
% 
% [X,Y,Z]=meshgrid(x,y,z);
% M=zeros(size(X));
% 
% c=1;
% for i=1:length(x)
%     for j=1:length(y)
%         for k=1:length(z)
%             M(i,j,k)=v(c);
%             c=c+1;
%         end
%     end
% end
%% new code for irregular grid
 % target regular grid for plotting
 xr = [257:2:383];
 yr = [257:2:383];
 zr = [1019:-2:1001];
 [Xr,Yr,Zr] = meshgrid(xr,yr,zr);
 
 Cr = zeros(size(Xr));
 
 % source data
 Cr = griddata(xt,yt,zt,ct,Xr,Yr,Zr) ;