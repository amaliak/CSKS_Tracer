
#include "bbfmm3d.hpp"
#include<iostream>
#include "mkl_vsl.h"

void SetMetaData(double& L, int& n, doft& dof, int& Ns, int& Nf, int& m, int& level, double& eps, int& nx, int& ny, int& nz) {
    

	ifstream fin;
	
	string filename = "./parameters.txt";

	fin.open(filename.c_str());
	
	if (!fin.good()){
		cerr << "Failed to open file " << filename << endl;
		throw runtime_error("Failed to open file!");
	}
    
	string line;
    getline(fin,line);
	getline(fin,line);
    
	line.erase(remove(line.begin(), line.end(), ' '),line.end());
    stringstream ss;
    ss << line;
    char comma;
    ss >> L >> comma >> n >> comma >> dof.s >> comma >> dof.f >> comma >> nx >> comma >> ny 
		>> comma >> nz >> comma >> m >> comma >> level >> comma >> eps;
    
	fin.close();

	Ns = nx*ny*nz;	// Number of sources in simulation cell
	Nf = nx*ny*nz;	// Number of field points in simulation cell
}

void read_xyz(const string& filenamex, int nx, const string& filenamey,int ny, const string& filenamez,int nz, vector3 *xyz) {
	ifstream fin;
	double* x = new double[nx];
	double* y = new double[ny];
	double* z = new double[nz];
	
	/* Read source */
	fin.open(filenamex.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenamex << endl;
		throw runtime_error("Failed to open file!");
	}
	for (int i=0;i<nx ;i++ )
	{
		fin >> x[i];
	}
	fin.close();

	fin.open(filenamey.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenamex << endl;
		throw runtime_error("Failed to open file!");
	}
	for (int i=0;i<ny ;i++ )
	{
		fin >> y[i];
	}
	fin.close();

	fin.open(filenamez.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filenamex << endl;
		throw runtime_error("Failed to open file!");
	}
	for (int i=0;i<nz ;i++ )
	{
		fin >> z[i];
	}
	fin.close();

	for (int k=0;k < nz;k++)
	{
		for (int j=0;j < ny ;j++ )
		{
			for (int i=0;i < nx ;i++)
			{
				xyz[k*nx*ny + j*nx + i].x = x[i];
				xyz[k*nx*ny + j*nx + i].y = y[j];
				xyz[k*nx*ny + j*nx + i].z = z[k];
			}
		}
	}


}

/*
 * Function: SetSources 
 * -------------------------------------------------------------------
 * read file from input
 */


void SetSources(vector3 *field, int Nf, vector3 *source, int Ns, double *q, int m,
                doft *dof, double L, int nx, int ny, int nz) {

	int l, i, j, k=0;

	/* Read source */
	string filename = "./input.bin";
	
	ifstream fin;
	
	fin.open(filename.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filename << endl;
		throw runtime_error("Failed to open file!");
	}

	fin.read((char*) q, m*Nf*dof->s*sizeof(double));
	fin.close();	

	read_xyz("./xcoord.txt",nx,"./ycoord.txt",ny,"./zcoord.txt",nz, source);

	for (i=0;i<Nf;i++) {
		field[i].x = source[i].x;
		field[i].y = source[i].y;
		field[i].z = source[i].z;
    }
}



int main(int argc, char *argv[]) {
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
    
    double L;       // Length of simulation cell (assumed to be a cube)
    int n;          // Number of Chebyshev nodes per dimension
    doft dof;
    int Ns;         // Number of sources in simulation cell
    int Nf;         // Number of field points in simulation cell
    int m;
    int level;
    double eps;
    int use_chebyshev = 1;
    int nx,ny,nz;

    SetMetaData(L, n, dof, Ns, Nf, m, level, eps, nx, ny, nz);

	vector3* source = new vector3[Ns]; // Position array for the source points
	vector3* field = new vector3[Nf];  // Position array for the field points

	double* q = new double[Ns*dof.s*m]; // Source array

	SetSources(field,Nf,source,Ns,q,m,&dof,L,nx,ny,nz);


	double err;
    double *stress      =  new double[Nf*dof.f*m];// Field array (BBFMM calculation)
	
	cout << "L    : " << L << endl;
	cout << "m    : " << Ns << endl;
	cout << "# chebyshev    : " << n << endl;
	cout << "eps            : " << eps << endl;
	cout << "level          : " << level << endl;
    cout << "# of cols in q : " << m << endl;

    
    /**********************************************************/
    /*                                                        */
    /*                 Fast matrix vector product             */
    /*                                                        */
    /**********************************************************/
    
    /*****      Pre Computation     ******/
    clock_t  t0 = clock();
    //cout << "Pre-computation start time: " << double(t0) / double(CLOCKS_PER_SEC) << endl;
    
	kernel_Exp Atree(&dof,L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree();
    
	clock_t t1 = clock();
    //cout << "Pre-computation end time: " << double(t1) / double(CLOCKS_PER_SEC) << endl;
    
	double tPre = t1 - t0;

    /*****      FMM Computation     *******/
    t0 = clock();
    cout << "FMM computation start time: " << double(t0) / double(CLOCKS_PER_SEC) << endl;
	
	H2_3D_Compute<kernel_Exp> compute(&Atree, field, source, Ns, Nf, q,m, stress);
    
	t1 = clock();
    double tFMM = t1 - t0;

	for (int i=0;i<Ns*dof.s*m ;i++ )
	{	
		q[i] = stress[i];
	}

	cout << "FMM2 computation starts " << endl;

    /*****      FMM Computation     *******/
    t0 = clock();
	
	H2_3D_Compute<kernel_Exp> compute2(&Atree, field, source, Ns, Nf, q,m, stress);
    
	t1 = clock();
    double tFMM2 = t1 - t0;

	for (int i=0;i<Ns*dof.s*m ;i++ )
	{	
		q[i] = stress[i];
	}
	
	cout << "FMM3 computation starts " << endl;

    /*****      FMM Computation     *******/
    t0 = clock();
	
	H2_3D_Compute<kernel_Exp> compute3(&Atree, field, source, Ns, Nf, q,m, stress);
    
	t1 = clock();
    double tFMM3 = t1 - t0;

    cout << "FMM computation end time: " << double(t1) / double(CLOCKS_PER_SEC) << endl;

    /*****      output result to binary file    ******/
    string outputfilename = "./stress.bin";
    write_Into_Binary_File(outputfilename, stress, m*Nf*dof.f);
    
    /**********************************************************/
    /*                                                        */
    /*              Exact matrix vector product               */
    /*                                                        */
    /**********************************************************/

	cout << "Pre-computation time: " << double(tPre) / double(CLOCKS_PER_SEC) << endl;
    cout << "FMM computing time:   " << double(tFMM) / double(CLOCKS_PER_SEC)  << endl;
	cout << "FMM total time:   "  << double(tPre+tFMM + tFMM2 + tFMM3) / double(CLOCKS_PER_SEC)  << endl;
    
    /*******            Clean Up        *******/
    
    delete []stress;
    return 0;
}
