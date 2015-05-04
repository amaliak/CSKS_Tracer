#
# File: Makefile
# Description: Makefile for BBFMM
# ------------------------------------------------------------------
# Use "make" to compile the code.
#
# Black-Box Fast Multipole Method (BBFMM)
# William Fong
# Stanford University
#
#

CC = icpc
LD = icpc

#CFLAGS  = -c -Wall -fpermissive -O3 -I ./include/ -I/usr/include -mavx -axAVX -xAVX -inline-forceinline -no-prec-div -ansi-alias  
#CFLAGS  = -c -O3 -mavx -axAVX -xAVX -no-prec-div -ansi-alias -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -I$(MKLROOT)/include/fftw -I ./include/ -I/usr/include 
#CFLAGS  = -c -O3 -mavx -axAVX -xAVX -no-prec-div -ansi-alias -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -I ./include/ -I/usr/include

CFLAGS  = -c -O3 -mavx -axAVX -xAVX -ipo -no-prec-div -ansi-alias -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -I ./include/ -I/usr/include

#LDPATH = -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -I$(MKLROOT)/include/fftw -L/ -L/usr/lib -I/usr/include 
LDPATH = -I $(MKLROOT)/include -L $(MKLROOT)/lib/intel64 -L/ -L/usr/lib -I/usr/include 

#LDFLAGS = -O3 -llapack -lblas -lrfftw -lfftw
LDFLAGS = -O3 -ipo -mkl -lrfftw -lfftw

MFLAGS = -lm 
#MFLAGS = $(MKLROOT)/lib/intel64/libfftw2xc_double_intel.a -lm 

PFLAG  =
SOURCES =  ./src/kernel_Types.cpp ./src/H2_3D_Tree.cpp ./src/read_metadata.cpp ./src/read_sources.cpp ./src/write_Into_Binary_File.cpp ./src/read_output.cpp

SOURCE = ./examples/get_input_through_routine_exp.cpp
SOURCE_FMM_EXP = ./examples/fmm_exp.cpp
SOURCE_FMM3_EXP = ./examples/fmm3_exp.cpp

OBJECT=$(SOURCES:.cpp=.o) $(SOURCE:.cpp=.o)
OBJECT_FMM_EXP=$(SOURCES:.cpp=.o) $(SOURCE_FMM_EXP:.cpp=.o)
OBJECT_FMM3_EXP=$(SOURCES:.cpp=.o) $(SOURCE_FMM3_EXP:.cpp=.o)

EXECUTABLE=  ./exec/get_input_through_routine_exp
EXECUTABLE_FMM_EXP=  ./exec/fmm_exp
EXECUTABLE_FMM3_EXP=  ./exec/fmm3_exp

all: $(SOURCES) $(SOURCE) $(EXECUTABLE)
$(EXECUTABLE): $(OBJECT)
	$(CC)  $(LDPATH) $(LDFLAGS)  $(OBJECT) -o $@ $(MFLAGS)

fmm_exp: $(SOURCES) $(SOURCE_FMM_EXP) $(EXECUTABLE_FMM_EXP)
$(EXECUTABLE_FMM_EXP): $(OBJECT_FMM_EXP)
	$(CC)  $(LDPATH) $(LDFLAGS)  $(OBJECT_FMM_EXP) -o $@ $(MFLAGS)

fmm3_exp: $(SOURCES) $(SOURCE_FMM3_EXP) $(EXECUTABLE_FMM3_EXP)
$(EXECUTABLE_FMM3_EXP): $(OBJECT_FMM3_EXP)
	$(CC)  $(LDPATH) $(LDFLAGS)  $(OBJECT_FMM3_EXP) -o $@ $(MFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o *~ ./src/*.o ./examples/*.o ./exec/fmm_exp ./exec/get_* ./output/*

tar:
	tar -zcvf BBFMM3D.tar.gz ./exec ./src ./include ./examples ./Makefile ./input ./output ./README.md
