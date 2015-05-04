class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int level, int n,  double epsilon, int
       use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
    virtual void setHomogen(string& kernelType,doft*dof) {
       homogen = 0;
       symmetry = 1;
       dof->f = 1;
       dof->s = 1;
       kernelType = "myKernel";}
       virtual void EvaluateKernel(vector3 fieldpos, vector3 sourcepos,
                               double *K, doft *dof) {
        double lx = 100; 
	double ly = 100;
	double lz = 100;
	double rx = (sourcepos.x - fieldpos.x)*(sourcepos.x - fieldpos.x)*(1.0/lx)*(1.0/lx);
	double ry = (sourcepos.y - fieldpos.y)*(sourcepos.y - fieldpos.y)*(1.0/ly)*(1.0/ly);
	double rz = (sourcepos.z - fieldpos.z)*(sourcepos.z - fieldpos.z)*(1.0/lz)*(1.0/lz);		
	double r =  sqrt( rx + ry + rz );
        double t0;         //implement your own kernel on the next line
       
	t0 = exp(-r*r);
       *K =  t0;
    }
};
