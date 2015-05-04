%% Compile mex code for new kernel
clear all
syms r;                          % distance of two points (radius)
%kernel = 1 ./ r.^2;               % radius basis function 
kernel = exp(-r^2);
outputfile = 'mexFMM3D';
homogen =  0;                    % K(ax, ay) = a^m K(x,y),=> homogen = m
symmetry = 1;                    % symmetric: 1; non-symmetric: 0; 
                                 % anti-symmetric: 0
make(r,kernel,homogen,symmetry,outputfile);