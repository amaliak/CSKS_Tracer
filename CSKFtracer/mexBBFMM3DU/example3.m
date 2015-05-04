% %% Compile mex code for new kernel
 clear all
 delete('./*.o');
 delete('./BBFMM3D/output/*.bin'); % make sure to delete output file if kernel is changed
 syms r;                           % distance of two points (radius)
% kernel = 1 ./ r.^2;               % radius basis function
 corlen = 10;
 kernel = exp(-(r.^2)./corlen);
% kernel = exp(-r.^2);
% kernel = exp(-abs(r));
 outputfile = 'mexFMM3D';
 homogen = 0;                    % K(ax, ay) = a^m K(x,y),=> homogen = m
 symmetry = 1;                    % symmetric: 1; non-symmetric: 0; 
%                                  % anti-symmetric: 0
 make(r,kernel,homogen,symmetry,outputfile);
%make(outputfile) 

%% Example: Q*H


% 3-D locations, stored column-wise i.e. (x | y | z)
%source = rand(Ns,3) - 0.5;
% field = rand(Nf,3) - 0.5; 
%field = source;

x = -62:4:62;
y = -62:4:62;
z = -9:3:9;

% x=0.001*x;
% y=0.001*y;
% z=0.001*z;
L = max(max(max(x) - min(x), max(y)-min(y)),max(z)-min(z));% Length of simulation cell (assumed to be a cube)

nx = length(x); ny = length(y); nz=length(z);
source = zeros(nx*ny*nz,3);

for i = 1:ny
    for j = 1:nx
        for k = 1:nz
            source((i-1)*nx*nz + (j-1)*nz + k ,1) = x(i);
            source((i-1)*nx*nz + (j-1)*nz + k ,2) = y(j);
            source((i-1)*nx*nz + (j-1)*nz + k ,3) = z(k);
        end
    end
end

%source = [0. 0. 0.;0. 0. 1.;0. 1. 0.;0. 1. 1.;1. 0. 0.;1. 0. 1.;1. 1. 0.;1. 1. 1.];
field = source;

% Info on dimensions
Ns  = nx*ny*nz;    % Number of sources in simulation cell
Nf  = Ns;    % Number of fields in simulation cell
m   = 1;        % number of columns of H

%H = randn(Ns,m); 
%H = ones(Ns,m); 
load H
% parameters
nCheb = 4;          % Number of Chebyshev nodes per dimension
level = 3;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation
%L =1;

%%
% Compute matrix-vectors product QH
% QH = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);
[QH,QH_exact] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);

%[QH] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);



 
 


