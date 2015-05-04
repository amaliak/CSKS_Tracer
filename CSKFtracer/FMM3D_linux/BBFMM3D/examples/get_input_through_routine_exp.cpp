
#include "bbfmm3d.hpp"
#include<iostream>
#include "mkl_vsl.h"

void SetMetaData(double& L, int& n, doft& dof, int& Ns, int& Nf, int& m, int& level, double& eps, int& nx, int& ny, int& nz) {
    L       = 91;    // Length of simulation cell (assumed to be a cube)
    n       = 6;    // Number of Chebyshev nodes per dimension
    dof.f   = 1;
    dof.s   = 1;

	nx      = 180;
	ny      = 180;
	nz      = 60;
	
	Ns      = nx*ny*nz;  // Number of sources in simulation cell
    Nf      = nx*ny*nz;  // Number of field points in simulation cell
    m       = 1;
    level   = 6;
    eps     = 1e-9;
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

	//cout << "x[0,1,..,nx-3:nx-1] = " << x[0] << ", "<< x[1] << ", " << x[nx-3] << ", " << x[nx-2] << ", " << x[nx-1] << endl;
	//cout << "y[0,1,..,nz-3:nz-1] = " << y[0] << ", "<< y[1] << ", " << y[ny-3] << ", " << y[ny-2] << ", " << y[ny-1] << endl;
	//cout << "z[0,1,..,nz-3:nz-1] = " << z[0] << ", "<< z[1] << ", " << z[nz-3] << ", " << z[nz-2] << ", " << z[nz-1] << endl;

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
 * Distributes an equal number of positive and negative charges uniformly
 * in the simulation cell ([-0.5*L,0.5*L])^3 and take these same locations
 * as the field points.
 */


void SetSources(vector3 *field, int Nf, vector3 *source, int Ns, double *q, int m,
                doft *dof, double L, int nx, int ny, int nz) {

	//double r[Ns]; /* buffer for random numbers */

	//VSLStreamStatePtr stream;
	//vslNewStream( &stream, VSL_BRNG_SFMT19937, 777 );
	//vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, Ns, r, 0.0, 1.0 );
    
	int l, i, j, k=0;
	

	/* Read source */
	
	string filename = "../input/randn.bin";
	
	ifstream fin;
	
	fin.open(filename.c_str(),ios::binary);
	if (!fin.good()){
		cerr << "Failed to open file " << filename << endl;
		throw runtime_error("Failed to open file!");
	}

	fin.read((char*) q, m*Nf*dof->s*sizeof(double));

	fin.close();	
	/*
	cout << "q[0:9]" << endl;
	
	for (int i=0;i<10 ; i++)
	{
		cout << q[i] << endl;
	}
	*/
	/*	
	// Distributes the sources randomly uniformly in a cubic cell
    for (l=0;l<m;l++) {
        for (i=0;i<Ns;i++) {
            for (j=0; j<dof->s; j++, k++){
				q[k] = 1.0;
            }
        }
    
    }
    */

	read_xyz("./xcoord.txt",nx,"./ycoord.txt",ny,"./zcoord.txt",nz, source);

//	for (i=1;i< 5 ;i++ )
//	{
//		cout << "i " << i << " x:" << source[i].x << " y:" << source[i].y << " z:" << source[i].z << endl;
//	}
//	
//	for (i=nx+1;i< nx+5 ;i++ )
//	{
//		cout << "i " << i << " x:" << source[i].x << " y:" << source[i].y << " z:" << source[i].z << endl;
//	}

    // Randomly set field points
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
    //vector3 source[Ns];    // Position array for the source points
    //vector3 field[Nf];     // Position array for the field points
    //cout << nx << " x " << ny << " x " << nz << endl;

	vector3* source = new vector3[Ns]; // Position array for the source points
	vector3* field = new vector3[Nf];  // Position array for the field points

    
	//double q[Ns*dof.s*m];  // Source array
    double* q = new double[Ns*dof.s*m]; // Source array

	SetSources(field,Nf,source,Ns,q,m,&dof,L,nx,ny,nz);


	double err;
    double *stress      =  new double[Nf*dof.f*m];// Field array (BBFMM calculation)
	
	cout << "L    : " << L << endl;
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
    cout << "Pre-computation start time: " << double(t0) / double(CLOCKS_PER_SEC) << endl;
    
	kernel_Exp Atree(&dof,L,level, n, eps, use_chebyshev);
    Atree.buildFMMTree();
    
	clock_t t1 = clock();
    cout << "Pre-computation end time: " << double(t1) / double(CLOCKS_PER_SEC) << endl;
    
	double tPre = t1 - t0;

    /*****      FMM Computation     *******/
    t0 = clock();
    cout << "FMM computation start time: " << double(t0) / double(CLOCKS_PER_SEC) << endl;
	
	H2_3D_Compute<kernel_Exp> compute(&Atree, field, source, Ns, Nf, q,m, stress);
    
	t1 = clock();
    double tFMM = t1 - t0;

    cout << "FMM computation end time: " << double(t1) / double(CLOCKS_PER_SEC) << endl;
    
	for (int i = 1;i < Ns ;i++ )
    {
		q[i] = stress[i];
    }
	
	/*****      FMM Computation     *******/
    t0 = clock();
    cout << "FMM computation2 start time: " << double(t0) / double(CLOCKS_PER_SEC) << endl;
	
	H2_3D_Compute<kernel_Exp> compute2(&Atree, field, source, Ns, Nf, q,m, stress);
    
	t1 = clock();
    double tFMM2 = t1 - t0;

    cout << "FMM computation2 end time: " << double(t1) / double(CLOCKS_PER_SEC) << endl;

    for (int i = 1;i < Ns ;i++ )
    {
		q[i] = stress[i];
    }

	/*****      FMM Computation     *******/
    t0 = clock();
    cout << "FMM computation3 start time: " << double(t0) / double(CLOCKS_PER_SEC) << endl;
	
	H2_3D_Compute<kernel_Exp> compute3(&Atree, field, source, Ns, Nf, q,m, stress);
    
	t1 = clock();
    double tFMM3 = t1 - t0;

    cout << "FMM computation3 end time: " << double(t1) / double(CLOCKS_PER_SEC) << endl;
	
	/*****   You can repeat this part with different source, field points and charges *****/
    
    /*vector3 source1[Ns];    // Position array for the source points
    vector3 field1[Nf];     // Position array for the field points
    double q1[Ns*dof.s*m];  // Source array
     
    SetSources(field1,Nf,source1,Ns,q1,m,&dof,L);
	*/     
    /****   Test interplation error   *****/
    
    //kernel_Exp testTree(&dof,1.0/4 ,2, n, eps, use_chebyshev);
	//double errtest = testInterplationErr(&testTree, 100, 100);
    
//	cout << "print out 5 vals from FMM" << endl;
//	cout.precision(15);
//	for (int i=0;i<10;i++) {
//		cout << double(stress[i]) << endl;
//    }

    
    /*****      output result to binary file    ******/
    string outputfilename = "../output/stress.bin";
    write_Into_Binary_File(outputfilename, stress, m*Nf*dof.f);
    
    /**********************************************************/
    /*                                                        */
    /*              Exact matrix vector product               */
    /*                                                        */
    /**********************************************************/
    /*
    // Field array (direct O(N^2) calculation)
	double *stress_dir  =  new double[Nf*dof.f*m];
	t0 = clock();
    DirectCalc3D(&Atree, field, Nf, source, q, m, Ns, &dof,0 , L, stress_dir);
	string outputfilename_dir = "../output/stress_dir.bin";
	write_Into_Binary_File(outputfilename_dir, stress_dir, m*Nf*dof.f);
    
	read_output(outputfilename_dir,stress_dir,Nf,m,dof);

	t1 = clock();
    double tExact = t1 - t0;
    cout << "print out 5 vals from direct calc" << endl;
    
	for (int i=0;i<5;i++) {
		cout << double(stress_dir[i]) << endl;
    }
	*/
	cout << "Pre-computation time: " << double(tPre) / double(CLOCKS_PER_SEC) << endl;
    cout << "FMM computing time:   " << double(tFMM) / double(CLOCKS_PER_SEC)  << endl;
    //cout << "FMM total time:   "  << double(tPre+tFMM) / double(CLOCKS_PER_SEC)  << endl;
    
	cout << "FMM total time:   "  << double(tPre+tFMM + tFMM2 + tFMM3) / double(CLOCKS_PER_SEC)  << endl;
    //cout << "Exact computing time: " << double(tExact) / double(CLOCKS_PER_SEC)  << endl;
    
    // Compute the 2-norm error
    //err = ComputeError(stress,stress_dir,Nf,&dof,m);
    //cout << "Relative Error: "  << err << endl;
    //err = ComputeError(stress,stress1,Nf,&dof,m);
    
    /*******            Clean Up        *******/
    
    delete []stress;
    //delete []stress_dir;
    //delete []stress1;
    return 0;
}
