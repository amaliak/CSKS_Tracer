%% Compile mex code for new kernel
clear all
% syms r;                          % distance of two points (radius)
% kernel = 1 ./ r.^2;               % radius basis function 
% %kernel = exp(-r^2);
outputfile = 'mexFMM3D';
% homogen =  -2;                    % K(ax, ay) = a^m K(x,y),=> homogen = m
% symmetry = 1;                    % symmetric: 1; non-symmetric: 0; 
%                                  % anti-symmetric: 0
make(outputfile);
 

%% Example: Q*H

% Info on dimensions
Ns  = 125;    % Number of sources in simulation cell
Nf  = 125;    % Number of fields in simulation cell
m   = 1;        % number of columns of H

% 3-D locations, stored column-wise i.e. (x | y | z)
%source = rand(Ns,3) - 0.5;
% field = rand(Nf,3) - 0.5; 
%field = source;

x = -0.5:0.25:0.5;
y = -0.5:0.25:0.5;
z = -0.5:0.25:0.5;

source = zeros(5*5*5,3);

for i = 1:5
    for j = 1:5
        for k = 1:5
            source((i-1)*25 + (j-1)*5 + k ,1) = x(i);
            source((i-1)*25 + (j-1)*5 + k ,2) = y(j);
            source((i-1)*25 + (j-1)*5 + k ,3) = z(k);
        end
    end
end

%source = [0. 0. 0.;0. 0. 1.;0. 1. 0.;0. 1. 1.;1. 0. 0.;1. 0. 1.;1. 1. 0.;1. 1. 1.];
field = source;
%H = rand(Ns,1); 
H = ones(Ns,1);

nCheb = 4;          % Number of Chebyshev nodes per dimension
L = 1.0;            % Length of simulation cell (assumed to be a cube)
level = 3;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation


% Compute matrix-vectors product QH
[QH,QH_exact] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);


 
 


