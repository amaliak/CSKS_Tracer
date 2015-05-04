/*!\file read_output.hpp
 read sources from binary files
*/
#ifndef __read_output_hpp__
#define __read_output_hpp__

#include"environment.hpp"
#include"bbfmm.h"
using namespace std;

void read_output(const string& filenameCharge, double *q, const int& Nf, const int& m, const doft& dof);

#endif //(__read_output_hpp__)
