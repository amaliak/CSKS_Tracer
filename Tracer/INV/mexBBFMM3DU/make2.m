function  make2(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Main Function%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% this file will call mexBBFMM3D.cpp
src1 = './BBFMM3D/src/H2_3D_Tree.cpp';
src2 = './BBFMM3D/src/kernel_Types.cpp';

disp(pwd)
eigenDIR = './eigen/';
fmmDIR = './BBFMM3D/include/';
mex('-O','./mexFMM3D.cpp',src1, src2,'-largeArrayDims',['-I',eigenDIR],['-I',fmmDIR],...
    '-llapack', '-lblas',...
    '-L/home/amaliak/fftw/lib', '-lfftw', '-lrfftw', '-lm','-g',...
    '-I/opt/intel/Compiler/11.1/084/Frameworks/mkl/include/fftw',...
    '-I/home/amaliak/fftw/include',...
    '-I.', '-output',filename)
disp('mex compiling is successful!')
end
