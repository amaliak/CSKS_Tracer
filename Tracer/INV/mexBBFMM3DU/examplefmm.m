function QH=examplefmm(H,nCheb)
% before using this function, use example3 and example3b to test FMM and create executables
% %% Compile mex code for new kernel
tic
clearvars -except H nCheb
%clear all
recompile = 0;
if recompile == 1
 %change kernelfun.hpp explicitly
 delete('./*.o');
 delete('./BBFMM3D/output/*.bin'); % make sure to delete output file if kernel is changed
 % syms r;                           % distance of two points (radius)
 % kernel = 1 ./ r.^2;               % radius basis function
 %corlen=70; 
 %kernel = exp(-(r^2)/corlen);
 outputfile = 'mexFMM3D';
 %homogen = 0;                    % K(ax, ay) = a^m K(x,y),=> homogen = m
 %symmetry = 1;                    % symmetric: 1; non-symmetric: 0; 
 %                                  % anti-symmetric: 0
 % make(r,kernel,homogen,symmetry,outputfile);
 % make2 file does not regenerate kernelfun.hpp, change it explicitly	
 make2(outputfile) 
end

%% Example: Q*H

load('./coord_htr.mat')
x_loc = x_htr - 382+(382-258)/2;
y_loc = y_htr - 382+(382-258)/2;
z_loc = z_htr - mean(z_htr);
lx = 60; lz=10;
z_loc = z_loc * (lx/lz);
L = max(max(max(x_loc) - min(x_loc), max(y_loc)-min(y_loc)),max(z_loc)-min(z_loc));% Length of simulation cell (assumed to be a cube)
source = [x_loc, y_loc,z_loc];
field = source;
% show cor lengths used as a reminder
fidl=fopen('./BBFMM3D/include/kernelfun.hpp');
for l=1:12
s=fgets(fidl);
end
disp('Displaying correlation lengths from kernelfun.hpp - Modify file directly if needed')
for l=1:3
s=fgets(fidl);
disp(s)
end
for l=1:6
s=fgets(fidl);
end
s=fgets(fidl);
disp('Covariance function')
disp(s)
fclose(fidl);

% Info on dimensions
Ns  = length(source);    % Number of sources in simulation cell
Nf  = Ns;    % Number of fields in simulation cell
m   = size(H,2);        % number of columns of H


if size(H,1)~=Ns
    disp('error in H')
end
%

% parameters
%nCheb = 4;          % Number of Chebyshev nodes per dimension
level = 5;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation
%L =1;
%%
% Compute matrix-vectors product QH
% QH = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);
%[QH,QH_exact] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);
QH = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);
toc

% exact


 
 

