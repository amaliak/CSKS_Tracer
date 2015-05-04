%%%%%%%%%%%%%% write ELEME and CONNE to MESH file%%%%%%%%%%%%%
clear all; close all;
% Geometry
m = 2025;
[x,y] = meshgrid(5:10:450,5:10:450);
x = x(:);y = y(:);
z = 5.*ones(size(y));

% Heterogeneous permeability field
load('CHl.mat');
rng(45)
px = CH'*randn(m,1); % permeability modifier is a Gaussian field
px = px - 1.1*min(px);
% plot permeability
px = reshape(px,45,45);
px = rot90(px); % rotate 90 degree counter clockwise
imagesc(px.*3e-12); colorbar; 
px = px(:); 


%% Modify ELEME file
% 1. write permeability modifier
% 2. impose Dirichilet boundary condition to certain cell

% non-BC cell
volx = 10*10*10; % volumn of each element
id = 46:1980;   % element ID is arranged in column-wise fashion
fid =fopen('ELEME', 'w');
fprintf(fid,'ELEME\n');
for i = 1:length(id)
    j = id(i);
    fprintf(fid,'%08d       ROCK1%10.4E%10.4E%10.4E%10.4E%10.4E%10.4E\n',j,volx,0,px(j),x(j),y(j),z(j));
end

% BC cell
volx = 10*10*10; % set the element as "inactive" by setting its volume to be zero
iBC = [1981:2025 1:45];
for i = 1:length(iBC)
j = iBC(i);
fprintf(fid,'%08d       ROCK1%10.4E%10.4E%10.4E%10.4E%10.4E%10.4E\n',j,volx,0,px(j),x(j),y(j),z(j));
end
fprintf(fid,'\n');

% 3. write CONNE
fprintf(fid,'CONNE\n');