function [U,S,V] = RandomizedCondSVDFMM(m,N)
tic
%Low-rank svd via randomized algorithm
%Performs fast SVD for covariance matrix and finds U, S, and V

% Arguments:
% m = size of covariance matrix number of unknowns
% N = rank of reduced rank svd 
% q is 1, or 2 (hardcoded below to = 1)
% Q is not in the input because it is calculated by BBFMM3D, see
% below
% set to 1 if you are changing covariance to reproduce mex file
flag_new_covariance = 1;
% To improve performance, randSVD is combined with BBFMM3D
% Care should be taken to use  BBFMM3D code consistent with the covariance
% matrix we are interested in - the Q matrix is calculated internally in
% the BBFMM3D code
%
% Amalia Kokkinaki
% modification of pkk code with in put from Harry J. Lee's
% function rnd_eig_fft.m and help from Ruoxi Wang on the use of mexBBFMM3D
% for mexBBFMM3D, the kernel information should be explicitly changed in file BBFMM3D/include/kernelfun.hpp
if N>m
    warning('N set to minimum of dimensions')
    N=m;
end
tic
% use "a" for randomized method; oversampling parameter
a = 9; 

%Y = (A*A')^q*A*Omega; % A(A(AO))) 3x(mlogm N)
% since A = Cov, A=A', Y=(Q^2q)*Q*Omega  
% for q = 1 Y = Q*Q*Q*Omega
q = 1;

%% setting up BBFMM3D
load('./coord_htr2.mat')
x_loc = x_htr2 - 382+(382-258)/2;
y_loc = y_htr2 - 382+(382-258)/2;
z_loc = z_htr2 - mean(z_htr2);
lx = 100 ; lz = 10;
z_loc = z_loc * (lx/lz);
L = max(max(max(x_loc) - min(x_loc), max(y_loc)-min(y_loc)),max(z_loc)-min(z_loc));
source = [x_loc, y_loc,z_loc];
field = source;

% Info on dimensions
Ns  = length(source);    % Number of sources in simulation cell
Nf  = Ns;    % Number of fields in simulation cell

% FMM parameters
nCheb = 4;          % Number of Chebyshev nodes per dimension
%level = 3;          % Level of FMM tree
level = ceil(log10(Ns));
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation


% create mex file for appropriate covariance function 

cd('mexBBFMM3DU/')
      
if flag_new_covariance == 1 
    % COMPILATION NEEDED IF KERNEL, LEVEL, OR NCHEB IS CHANGED
    disp('New covariance kernel being created')
    % NOTE: IF TYPE OF KERNEL IS CHANGED, IT NEEDS TO BE CHANGED IN KERNELFUN.HPP
    %check = input('Did you change kernelfun.hpp');
    %if ~check
    %    disp('Change kernelfun.hpp now and hit space');
    %    pause;
    %else
    %    disp('Covariance info is OK, double-check with print out below.')
    %end
    % see examplefmm.m and example3b.m for details on mex compilation
    delete('./*.o');
    delete('./BBFMM3D/output/*.bin'); % make sure to delete output file if kernel is changed
    outputfile = 'mexFMM3D';
    make2(outputfile)
else
    disp('Using covariance kernel from last compilation')
    disp('Double-check covariance info is OK with print-out below')
end

% printout for covariance info
% show cor lengths used as a reminder
fidl=fopen('./BBFMM3D/include/kernelfun.hpp');
for l=1:12; s=fgets(fidl); end;
disp('Displaying correlation lengths from kernelfun.hpp - Modify file directly if needed')
for l=1:3; s=fgets(fidl); disp(s); end;
for l=1:6; s=fgets(fidl); end;
s=fgets(fidl);
disp('Covariance function')
disp(s)
fclose(fidl);
%pause; 

% example usage [QH] = mexFMM3D(source, field,H,nCheb, level, L, use_chebyshev);
% run a single mat-vec multiplicaiton to see the accuracy of FMM
disp('Running test for accuracy')
disp('Also creating the Pre-computation files')
%[QH QHexact]=mexFMM3D(source, field,Omega(:,1),nCheb, level, L, use_chebyshev);
testOmega=randn(m,1);
QH = mexFMM3D(source,field,testOmega(:,1),nCheb,level,L,use_chebyshev);
disp('FMM error is not being computed')
check = 1; %input('Is accuracy good enough?')
while ~check
	disp(['Current number of Cheb nodes is',num2str(nCheb)])
	nCheb = input('New number of Cheb nodes is');
	disp('Running test for accuracy')
	% recompiling..
	delete('./*.o');
        delete('./BBFMM3D/output/*.bin'); % make sure to delete output file if kernel is changed
        outputfile = 'mexFMM3D';
        make2(outputfile)
	[QH QHexact]=mexFMM3D(source, field,testOmega(:,1),nCheb, level, L, use_chebyshev);
	check = input('Is accuracy good enough?')
end	
clear  QH QHexact

%% Randomized SVD 

disp('Starting randSVD')
% the following is a multiplication of an mxm matrix with an thin matrix
% FMM is used to make it faster
% note that Q is provided by mexFMM3D, for each different Q need to
% recompile the mex file

% find how many processors are available
disp('Setting up parfor')
% chekck if there is existing parallel job running
% if so, delete it
p = gcp('nocreate');
if ~isempty(p)
    delete(gcp)
end

myCluster  = parcluster('local');
delete(myCluster.Jobs);

% start parallel job, find out number of workers available
poolobj = parpool('local');
noproc = poolobj.NumWorkers;
%noproc=12;
subH = (N+a)/noproc;
if mod(N+a,noproc)~=0
    disp(['N is ',num2str(N)])
    disp(['a is ',num2str(a)])
    disp(['noproc is ',num2str(noproc)])
    disp('Adjusting a in randSVD to fit number of processors')
    %decrease a to fit noproc
    % no need to regenerate Omega, just use fewer columns
    a = floor(subH)*noproc - N; 
    if a < 1
	a = ceil(subH)*noproc-N;
    end
    if (N+a) < (N+1)
	display('Error in randomized svd');
	keyboard;
    end
    disp(['New a is ',num2str(a),' and N+a is ',num2str(N+a)])		
    %a = input('Adjust ling a to fit noproc');
end
      
subH = (N+a)/noproc;

if mod(subH,1) > 0; keyboard; end;
      
if ~exist('Omega.mat')
      Omega = randn(m,N+a);
      save Omega.mat Omega
else
      load Omega.mat Omega
end
      
      
Y=zeros(m,N+a);    

% Y = Q*Omega;
parfor i=1:noproc
    start=subH*(i-1)+1;
    fin=i*subH;
    Hi=Omega(:,start:fin);
    res(:,:,i)=mexFMM3D(source, field,Hi,nCheb, level, L, use_chebyshev);
    disp(['Multiplication of colums from ',num2str(start),' to ',num2str(fin)])
end
disp('Computation of Y, Step 1/5 finished');

for i=1:noproc
    start=subH*(i-1)+1;
    fin=i*subH;
    Y(:,start:fin)=res(:,:,i);
end
% Y = (Q*Q)^q*(Q*Omega);
for d = 1 : 2*q
   parfor i=1:noproc
    	start=subH*(i-1)+1;
    	fin=i*subH;
   	Hi=Y(:,start:fin);
    	res(:,:,i)=mexFMM3D(source, field,Hi,nCheb, level, L, use_chebyshev);
    	disp(['Multiplication of colums from ',num2str(start),' to ',num2str(fin)])    
    end
    disp(['Computation of Y, Step ',num2str(d+1),'/5 finished'])
    
    for i=1:noproc
    	 start=subH*(i-1)+1;
     	fin=i*subH;
     	Y(:,start:fin)=res(:,:,i);
    end    
end
% step 2
[R,~,~] =svd(Y,0); 

% step 3  B = R'*Q --> B' = Q'R = QR
disp(['Computation of B, Last Step Started']);
parfor i=1:noproc
        start=subH*(i-1)+1;
        fin=i*subH;
        Hi=R(:,start:fin);
        res(:,:,i)=mexFMM3D(source, field,Hi,nCheb, level, L, use_chebyshev);
	disp(['Multiplication of colums from ',num2str(start),' to ',num2str(fin)])
end

disp(['Computation of B, Final Step finished']);
for i=1:noproc
    start=subH*(i-1)+1;
    fin=i*subH;
    B(:,start:fin)=res(:,:,i);
end

% step 4
[V,S,Ut] = svd(B,0);

% step 5
U = R*Ut;
U = U(:,1:N); S = S(1:N,1:N); V = V(:,1:N);

%TestError = norm(U*S*V'-A)/norm(A)

p = gcp('nocreate');
if ~isempty(p)
    delete(gcp)
end

cd ..

toc
end
