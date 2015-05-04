/*!\file	read_output.cpp
 source file to field, source and charge information from binary files.
*/

#include"read_output.hpp"
#include"environment.hpp"
using namespace std;

void read_output(const string& filenameCharge, double *q, const int& Nf, const int& m, const doft& dof) {
    ifstream fin;

    /* Read Output */
    fin.open(filenameCharge.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenameCharge << endl;
		throw runtime_error("Failed to open file!");
	}
    fin.read((char*) q, m*Nf*dof.f*sizeof(double));
    fin.close();
}