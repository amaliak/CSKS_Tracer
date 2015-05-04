function prod = callfmm(thinmat,flag1)

% function that performs the multiplication of a covariance
% with a thin matrix thinmat with the 3D FMM method (matlab implementation)
% the covariance is implemented internally within the FMM 
% make sure that the parameters of the covariance are consistent withi those in the
% calling function
% Note that in the original usage, because fmm was used before, the covariance parameters in the
% c++ file are from there and should be ok

% par for is utilized to split the job to available processors or cores

% prod = fatmat * thinmat
% flag1 --> if true, before the multiplication a comparison is performed with the direct method
	
% Setting up parameters
m = size(thinmat,1);
N = size(thinmat,2);

% parameters for covariance
% coordinates
load('./coord_htr2.mat')
x_loc = x_htr2 - 382+(382-258)/2;
y_loc = y_htr2 - 383+(382-258)/2;
z_loc = z_htr2 - mean(z_htr2);
% scaling for anisotropy
lx = 100; lz = 10; 
z_loc = z_loc * (lx/lz);
L=max(max(max(x_loc)-min(x_loc),max(y_loc)-min(y_loc)),max(z_loc)-min(z_loc));
source = [x_loc, y_loc, z_loc];
field = source; 
Ns = m; 
Nf = m; 
nCheb = 4; 
level = ceil(log10(Ns));
use_chebyshev = 1; 

cd('mexBBFMM3DU/')
% compilation of mex file for generation of precomputation files
if flag1
	disp('Recompiling and testing for accuracy')
	delete('./*.o');
	delete('./BBFMM3D/output/*.bin');
	outputfile = 'mexFMM3D';
	make2(outputfile)
	[QH QHexact]=mexFMM3D(source, field, thinmat(:,1),nCheb,level,L,use_chebyshev);
	check = input('Is accuracy good enough? ');
	while ~check
		disp(['Current number of Cheb nodes is',num2str(nCheb)])
		nCheb = input('New number of Cheb nodes is  ');
		disp('Recompiling and Running test for accuracy');
		delete('./*.o');
		delete('./BBFMM3D/output/*.bin');
		outputfile = 'mexFMM3D';
		make2(outputfile)
		[QH QHexact]=mexFMM3D(source,field, thinmat(:,1),nCheb,level,L,use_chebyshev);
		check = input('Is accuracy good enough? ');
	end
	clear QH QHexact
else
	disp('Recompiling')
	delete('./*.o');
	delete('./BBFMM3D/output/*.bin');
	outputfile = 'mexFMM3D';
	make2(outputfile)
	QH = mexFMM3D(source,field,thinmat(:,1),nCheb,level,L,use_chebyshev);
	clear QH
end

% parallelize actual multiplications
% Set up parfor
disp('Setting up parfor')
p = gcp('nocreate');
if ~isempty(p)
	delete(gcp)
end

myCluster = parcluster('local');
delete(myCluster.Jobs);
poolobj = parpool('local');
noproc = poolobj.NumWorkers;
if N < noproc
	noproc = N;
end 
setenv('LD_LIBRARY_PATH','/home/amaliak/fftw/lib')
subH = N / noproc; 
% start gives the starting index for each par for loop
% if start = [1 5], we have two(=size(start)) parfor loops
% one for columns 1:4 and one for columns 5:end 
start = [1:subH:N];
start = [floor(start) N+1]; %gives indices for splitting the columns
if length(start)~=noproc+1; disp('error'); keyboard; end; 
parfor i = 1:noproc
	col_1=start(i);
	col_n=start(i+1)-1;
	H = thinmat(:,col_1:col_n);
	res{i}=mexFMM3D(source,field,H,nCheb,level,L,use_chebyshev);
	disp(['Multiplication of columns from ',num2str(col_1),' to ',num2str(col_n)])
end

for i=1:noproc
	prod(:,start(i):start(i+1)-1) = res{i};
end 

disp('Matrix vector multiplication done')

% testing
%test = mexFMM3D(source,field,thinmat,nCheb,level,L,use_chebyshev);
cd ..
