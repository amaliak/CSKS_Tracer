% %% Compile mex code for new kernel
 clear all
 delete('./*.o');
 delete('./BBFMM3D/output/*.bin'); % make sure to delete output file if kernel is changed
 syms r;                           % distance of two points (radius)
% kernel = 1 ./ r.^2;               % radius basis function 
 corlen=80;
 kernel = exp(-(r^2)/corlen); 
 outputfile = 'mexFMM3D';
 homogen = 0;                    % K(ax, ay) = a^m K(x,y),=> homogen = m
 symmetry = 1;                    % symmetric: 1; non-symmetric: 0; 
%                                  % anti-symmetric: 0
 make(r,kernel,homogen,symmetry,outputfile);
%make(outputfile) 

%% Example: Q*H

% Info on dimensions
Ns  = 9261;    % Number of sources in simulation cell
Nf  = 9261;    % Number of fields in simulation cell
m   = 1;        % number of columns of H

% 3-D locations, stored column-wise i.e. (x | y | z)
%source = rand(Ns,3) - 0.5;
% field = rand(Nf,3) - 0.5; 
%field = source;

x = -0.5:0.05:0.5;
y = -0.5:0.05:0.5;
z = -0.5:0.05:0.5;
x=100*x;
y=100*y;
z=100*z;

source = zeros(21*21*21,3);

for i = 1:21
    for j = 1:21
        for k = 1:21
            source((i-1)*21*21 + (j-1)*21 + k ,1) = x(i);
            source((i-1)*21*21 + (j-1)*21 + k ,2) = y(j);
            source((i-1)*21*21 + (j-1)*21 + k ,3) = z(k);
        end
    end
end

%source = [0. 0. 0.;0. 0. 1.;0. 1. 0.;0. 1. 1.;1. 0. 0.;1. 0. 1.;1. 1. 0.;1. 1. 1.];
field = source;

%H = ones(Ns,m); 
load H
% parameters
nCheb = 4;          % Number of Chebyshev nodes per dimension
level = 3;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation
L =100;

%%
% Compute matrix-vectors product QH
% QH = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);
[QH,QH_exact] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);

%[QH] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);



 
 


