function make(filename)
% Produce mex file with 'filename' for the input kernel. File extension is automatically added. 
% For each new kernel type, run this make file first before call mexFMM2D. 
% The new kernel will be embedded in a new file 'kernelfun.hpp'

% Usage:
% syms r; f = exp(-abs(r)/30); filename = 'expfun';
% make(r,f,filename)

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Main Function%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% convert kernel to C readable format
% var = char(var);
% 
% ckernel = ccode(kernel);
% 
% % Pass the Kernel to kernelfun.hpp
% fid = fopen('./BBFMM3D/include/kernelfun.hpp','w+');
% fprintf(fid,'class myKernel: public H2_3D_Tree {\n');
% fprintf(fid,'public:\n');
% fprintf(fid,'    myKernel(doft* dof, double L, int level, int n,  double epsilon, int\n');
% fprintf(fid,'       use_chebyshev):H2_3D_Tree(dof,L,level,n, epsilon, use_chebyshev){};\n');
% fprintf(fid,'    virtual void setHomogen(string& kernelType) {\n');
% fprintf(fid,'       homogen = %d;\n', homogen);
% fprintf(fid,'       symmetry = %d;\n', symmetry);
% fprintf(fid,'       kernelType = "myKernel";}\n');
% fprintf(fid,'    virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,\n');
% fprintf(fid,'                               double *K, doft *dof) {\n');
% fprintf(fid,'        double %s =  sqrt((sourcepos.x - fieldpos.x)*(sourcepos.x - fieldpos.x) + (sourcepos.y - fieldpos.y)*(sourcepos.y - fieldpos.y) + (sourcepos.z - fieldpos.z)*(sourcepos.z - fieldpos.z));\n',var);
% fprintf(fid,'        double t0;         //implement your own kernel on the next line\n');
% fprintf(fid,'       %s\n',ckernel);
% fprintf(fid,'       *K =  t0;\n');
% fprintf(fid,'    }\n');
% fprintf(fid,'};\n');
% fclose(fid);

% this file will call mexBBFMM3D.cpp
src1 = './BBFMM3D/src/H2_3D_Tree.cpp';
src2 = './BBFMM3D/src/kernel_Types.cpp';

disp(pwd)
eigenDIR = './eigen/';
incDIR1 = '/usr/include';
fmmDIR = './BBFMM3D/include/';
mex('-v','./mexFMM3D.cpp',src1, src2,'-largeArrayDims',...
    ['-I',eigenDIR],['-I',fmmDIR],['-I',incDIR1],...
    '-L/usr/lib','-I/usr/include',...
    '-llapack', '-lblas','-L/home/amaliak/fftw/lib','-lrfftw', '-lfftw',...
    '-I/opt/intel/Compiler/11.1/084/Frameworks/mkl/include/fftw',...
    '-I/home/amaliak/fftw/include',...
    '-I.','-output',filename)
disp('mex compiling is successful!')
end


%'-I/opt/intel/composer_xe_2013.5.192/mkl/include/fftw',...
