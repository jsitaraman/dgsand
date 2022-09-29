/* Discontinuous Galerkin Sand Box for testing algorithms
 * efficiency and other good stuff. 
 * Jay Sitaraman 04/05/2021
 *
 * compile as 
 * gcc dgsand.c -lm -o dgsand
 *
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
extern "C" {
  void parseInputs(char *inputfile,char *gridfile,int *pde, int *itype, int *nsteps, double *dt, int *nsave,
		   int *ireg);
  void output_params();
  void readgrid2D(char *gridfile, double **xcoord, int **elem2node,int **ibc,
		  int *p,  int *nnodes, int *nelem, int *nbnodes);
  void find_faces(int *bface,
		  int **elem2face,
		  int **faces,
		  int *nfaces,
		  int *ibc,
		  int nsurfcells,
		  int nbnodes,
		  int nv, 
		  int nvmax); 
  
  void init_data(double *x,double *q, double *X, double *Q, int d, int nfields, int e, int p);  
  void INIT_FIELDS(double *xcoord,int *e2n, 
		   double *Q, double *x, double *q, // data populated here
		   int *iptr, int pde, int e, int p,
		   int d, int nbasis, int itype, int nelem, int pc,int imesh);  
  int number_of_fields(int pde,int d);
  int order_to_basis(int etype, int p);
  int get_ngElem(int etype, int p);
  int get_ngGL(int etype, int p);
  int face_per_elem(int etype);
  
  int FIND_NECUT(double x0, double *x,int* iptr, int d, int e, int p, int nelem, int pc, int imesh);

  
  double total_area(double *detJ, int etype, int p, int d, int nelem);
  void CUT_CELLS(double x0, double *x, double* xcut, int* iptr, int* cut2e, int d, int e, int p,
		 int nelem, int pc, int *cut2face, int* cut2neigh, int* elem2face, int* faces,
		 int* iblank, int* cutoverset, int imesh, int ng);
  void MASS_MATRIX(double *mass,double *x, int *iptr, int d, int e, int p, int nelem, int pc, int imesh);
  void CUT_MASS_MATRIX(double *mass,double *x, double *Jinv, int *iptr, double *xcut,
		       double *detJcut, int *iptrc, int d, int e, int p, int nelem,
		       int pc, int pccut, int necut, int* cut2e, int imesh);

  void COMPUTE_GRID_METRICS(double *x, double *bv, double *bvd,double *JinvV, 
			    double *detJ,double *bf, double *bfd, double *JinvF, double *faceWeight,
			    int *iptr, int d, int e, int p, int nelem, int pc, int* iblank, int imesh);

  void COMPUTE_CUT_METRICS(double *x, double *JinvV, 
			   double *detJ,double *JinvF,
			   int *iptr, int d, int e, int p, int pc, int pccut,
			   double *xcut, double *bvcut, double *bvdcut, 
			   double *JinvVcut,double *JinvFcut,
			   double *detJcut, double *bfcutL, double *bfdcutL,  
			   double *bfcutR, double *bfdcutR,  
			   double *fwcut, int* iptrc, int necut, int* cut2e,
			   int* cut2neigh, int imesh, int* cutoverset);
  
  void COMPUTE_RHS(double *R,double *mass,double *bv, double *bvd, double *JinvV, double *detJ,
		 double *bf, double *bfd, double *JinvF,
		 double *faceWeight, double *fnorm, double *fflux,
		 double *x, double *q, int *elem2face, int *iptr, int *iptrf, int *faces,
		 int pc, int pf, int pccut, int pde, int d , int e, int p, int nfaces, int nelem,
                 double *bvcut, double *bvdcut,double *detJcut,
                 double *bfcutL, double *bfcutR,double *fwcut, double *fcflux,
                 int* OSFnseg, int* OSFeID, double* OSFxn, double* OSFshpL, double* OSFshpR, double* OSFflux,
		 int *iptrc, int necut, int* cut2e, int *cut2face, int* cut2neigh, int* iblank, int ireg,
		 int *cutoverset,int imesh);
  void UPDATE_DOFS(double *qdest, double coef, double *qsrc, double *R, int ndof);

  void OUTPUT_TECPLOT(int meshid, int step,double *x, double *q,
		    int pc, int *iptr, int pde, int d, int e, int p, int nelem);

  void SETUP_OVERSET(int* cut2e, int* cut2eB, int* cut2neighA, 
		     int* iptrA, int* iptrB, int* iptrcA, int* iptrcB, 
		     double* xA, double* xB, double* xcutA, double* xcutB,
                     double* bfcutLA, double* bfcutRA, double* JinvA, double* JinvB,
		     int* OSFnseg, int* OSFeID, double* OSFxn, double* OSFshpL, double* OSFshpR,
		     int d, int e, int p, int pc, int pccut, int necutA, int necutB, int nelemB); 

  void EXCHANGE_OVERSET(double* OSFflux, double* OSFshpL, double* OSFshpR, 
                        int* OSFnseg, int* OSFeID, int* cut2eA,
                        double* qA, double* qB, int* cutoverset, int* iptrcA, 
                        int* iptrA, int* iptrB, int necutA, int pccut, 
                        int d, int e, int p, int pc, int pde, int imesh);

}

#include<vector>

class dgsand
{
 private:

  /// Mass matrix
  std::vector<double> mass;
  /// face normals and face flux
  std::vector<double> fnorm,fflux;
  /// face cut flux
  std::vector<double> fcflux;
  /// pde type, only Euler eqns are implemented now
  int pde;
  /// number of fields of this pde
  int nfields;
  /// initialization type
  int itype;
  ///  use Tikhonov reg. on cut cells, 0: don't
  int ireg; 			    
  /// overset inputs
  int nmesh = 2;
  
  /// Mesh geometry inputs
  /// number of nodes in grid
  int nnodes;
  /// number of primal nodes in physical boundaries
  int nbnodes;
  /// number of faces
  int nfaces;
  /// boundary condition node indices and their types
  int *ibc;
  /// faces connectivity graph
  int *faces;
  /// coordinates of the original grid
  double *xcoord;
  /// element to node connectivity
  int *elem2node;
  /// element to face connectivity
  int *elem2face;
  
  /* class local variables */
  int pf,ndof,nsteps,nsave;
  double dt;

 public:
  /// field quantities
  std::vector<double> q,qstar,Q,R;
  /// volumetric basis quantities and grid
  std::vector<double> x,bv,bvd,JinvV,detJ;
  /// face basis quantities
  std::vector<double> bf,bfd,JinvF,faceWeight;
  /// number of elements
  int nelem;
  /// pointer in to data arrays
  std::vector<int> iptr,iptf;
  int pc; 
  /// number of bases for solution
  int nbasis;
  /// number of bases for grid
  int nbasisx;

  /// Cut quantities
  /// cut face information
  std::vector<double> xcut,bvcut,bvdcut,JinvVcut,detJcut;
  /// cut face quantities
  std::vector<double> bfcutL,bfdcutL,JinvFcut,fwcut,bfcutR,bfdcutR;
  /// cutface connecitivity and iblanking
  std::vector<int> cut2e,iblank, cutoverset;  
  /// pointer in to data arrays for cut cells
  std::vector<int> iptrc;
  int pccut;
  /// cut face connectivity and neighbors
  std::vector<int> cut2face,cut2neigh;    
  /// number of edges cut
  int necut;
  /// Overset Face (OSF): number of segments for each cut cell and 
  /// neighbor cell ID for each segment gauss point
  std::vector<int> OSFnseg, OSFeID;
  /// Overset Face (OSF): shape function values at each segment gauss pt
  std::vector<double> OSFshpL, OSFshpR, OSFxn,OSFflux;
  
  /// dimensions (only implemented for 2D now)
  const int d=2;
  /// p order
  int p;
  /// element type
  int etype=0;   

  /// generic constructor
  dgsand(): p(1),ndof(0),nfields(4),pde(1),necut(0),etype(0){};
  /// generic destructor
  ~dgsand() { free(xcoord); free(elem2node); free(ibc); free(elem2face); free(faces);}
  /// setup the case using input file
  void setup(char *input_file)
    {
        /* parse inputs */
	char grid_file[20];
	parseInputs(input_file,grid_file, &pde,&itype,&nsteps,&dt,&nsave,&ireg);
	nfields=number_of_fields(pde,d);
	
	/* read a 2D grid */
	readgrid2D(grid_file,&xcoord,&elem2node,&ibc,&p,&nnodes,&nelem,&nbnodes);
	nbasis=order_to_basis(etype,p);        // basis for solution
	nbasisx=order_to_basis(etype,p+(p==0));// basis for grid
	
	/* find element to face connectivity */
	find_faces(elem2node,&elem2face,&faces,&nfaces,ibc,nelem,nbnodes,3,nbasisx);
	
	/* allocate memory */
	/* field parameters per element */
	int qsize=(nfields*nbasis*nelem);
	q.resize(qsize);     // modal coefficients		  
	qstar.resize(qsize); // modal coefficients		  
	Q.resize(qsize);	   // values at physical locations  
	R.resize(qsize);	   //solution residual             
	
	/* geometrical parameters per volume QP of each element */
	/* TODO: some of these such as bv and JinvV can be optimized or omitted */
	int ngElem=get_ngElem(etype,p);
	int ngGL=get_ngGL(etype,p);
	
	x.resize((d*(nbasisx))*nelem);                 // coord modal coefficients (p0 mod)   
	bv.resize((nbasis*ngElem));	       // basis value at volume QP	     
	bvd.resize((d*nbasis*ngElem*nelem)); // basis derivative value at volume QP 
	JinvV.resize((d*d*ngElem*nelem));    // J^{-1} at volume QP		     
	detJ.resize((ngElem*nelem));	       // |J| at volume QP		     
	iblank.resize(nelem);			       // iblank array                    
	/* geometrical parameters per face QP of each element */
	/* TODO : some these such as bf and JinvF can optimized/omitted */
	int fpe = face_per_elem(etype);
	bf.resize(nbasis*ngGL*fpe);            // basis value at face QP	  
	bfd.resize(d*nbasis*ngGL*fpe*nelem);   // basis der. value at face QP  
	JinvF.resize(d*d*ngGL*fpe*nelem);      // J^{-1} at face QP		  
	faceWeight.resize(d*ngGL*fpe*nelem);   // faceNormals at face QP	  
	mass.resize(nbasis*nbasis*nelem);		       // mass matrix      
	fnorm.resize(d*ngGL*nfaces);	       // face normals		                                  
	fflux.resize(3*nfields*ngGL*nfaces);   // face fields and flux               
	/* pointer array into each data array above */
	pc=11; // number of unique sizes with elements
	pf=2;  // number of unique sizes associated with faces
	iptr.resize(pc*nelem);
	iptf.resize(pf*nfaces);      
	/* set the pointers, TODO: this has to change when there is a variety of elements */
	for(int i=0;i<nelem;i++){    
	  int ix=pc*i;
	  iptr[ix]+=i*(nfields*nbasis);               // q, Q, R
	  iptr[ix+1]+=i*(d*(nbasisx));                // x
	  iptr[ix+2]+=0;                              // bv (this is same per element type)
	  iptr[ix+3]+=i*(d*nbasis*ngElem);  // bvd
	  iptr[ix+4]+=i*(d*d*ngElem);       // JinvV
	  iptr[ix+5]+=i*(ngElem);           // detJ
	  
	  iptr[ix+6]+=0;                              // bf (this is same per element type)
	  iptr[ix+7]+=i*(d*nbasis*ngGL*fpe);// bfd
	  iptr[ix+8]+=i*(d*d*ngGL*fpe);     // JinvF
	  iptr[ix+9]+=i*(d*ngGL*fpe);       // faceWeight
	  iptr[ix+10]+=i*(nbasis*nbasis);             // mass 
	}
	for(int i=0;i<nfaces;i++){
	  int ix=pf*i;
	  iptf[ix]+=(i*d*ngGL);            //faceNormal
	  iptf[ix+1]+=(i*3*nfields*ngGL);  //faceFlux
	}
      };

    void init(int imesh) {

      /* initialize fields on all the elements */
      INIT_FIELDS(xcoord,
		  elem2node,
		  Q.data(),
		  x.data(),
		  q.data(),
		  iptr.data(),
		  pde,etype,p,d,nbasis,itype,nelem,pc,imesh);
      
      /* compute grid metrics */
      COMPUTE_GRID_METRICS(x.data(),
			   bv.data(),
			   bvd.data(),
			   JinvV.data(),
			   detJ.data(),
			   bf.data(),
			   bfd.data(),
			   JinvF.data(),
			   faceWeight.data(),
			   iptr.data(),
			   d,etype,p,nelem,pc,iblank.data(),imesh);      
    }

    void cut(double x0,int imesh)
    {
      int ngElem=get_ngElem(etype,p);
      int ngGL=get_ngGL(etype,p);
      int maxseg = 20; 
      necut = FIND_NECUT(x0,x.data(),iptr.data(),d,etype,p,nelem,pc,imesh);
      if (necut > 0)  {
	//create all the cut cell pointers
	xcut.resize(d*3*necut);
	cut2e.resize(necut);       
	cut2face.resize(necut*3);  // map between cut face and orig face id                       
	cut2neigh.resize(necut*3); // map between cut face and R side neighbor
	cutoverset.resize(necut*3);        // cutoverset array 

	bvcut.resize(necut*nbasis*ngElem);         // basis value at volume QP	   
	bvdcut.resize(necut*d*nbasis*ngElem);// basis derivative value at volume QP
	JinvVcut.resize(necut*d*d*ngElem);   // J^{-1} at volume QP		   
	JinvFcut.resize(necut*d*d*ngElem);   // J^{-1} at volume QP		   
	detJcut.resize(necut*ngElem);	       // |J| at volume QP                   

	int fpe=face_per_elem(etype);
	bfcutL.resize(nbasis*ngGL*fpe*necut);       // basis value at face QP      
	bfdcutL.resize(d*nbasis*ngGL*fpe*necut);    // basis der. value at face QP 
	bfcutR.resize(nbasis*ngGL*fpe*necut);	      // basis value at face QP      
	bfdcutR.resize(d*nbasis*ngGL*fpe*necut);    // basis der. value at face QP 
	fwcut.resize(d*ngGL*fpe*necut);	      // faceNormals at face QP      
	fcflux.resize(3*nfields*ngGL*fpe*necut);    // face fields and flux        

	// overset face quantities
        OSFnseg.resize(necut);
        OSFeID.resize(necut*maxseg*ngGL); 		// neigh mesh cell ID for each overset QP
        OSFxn.resize(necut*maxseg*ngGL*d); 		// L*norm for each overset QP
        OSFshpL.resize(necut*maxseg*ngGL*nbasis); 	// neigh mesh basis value at each overset QP
        OSFshpR.resize(necut*maxseg*ngGL*nbasis); 	// neigh mesh basis value at each overset QP
        OSFflux.resize(necut*maxseg*ngGL*3*nfields);	// flux container for each overset QP
printf("OSFflux size = %i\n",necut*maxseg*ngGL*3*nfields); 	
	pccut = 17; 
	printf("nbasisx = %i\n",nbasisx); 
	iptrc.resize(pccut*necut);
	
	for(int i=0;i<necut;i++){
	  int ix=pccut*i;
	  iptrc[ix]+=i*(nfields*nbasis);               // q, Q, R
	  iptrc[ix+1]+=i*(d*3);                // x
	  iptrc[ix+2]+=i*(nbasis*ngElem);   // bv (this is NOT same per element type)
	  iptrc[ix+3]+=i*(d*nbasis*ngElem);  // bvd
	  iptrc[ix+4]+=i*(d*d*ngElem);       // JinvV
	  iptrc[ix+5]+=i*(ngElem);           // detJ
	  
	  iptrc[ix+6]+=i*(nbasis*ngGL*fpe); // bf (this is NOT same per element type)
	  iptrc[ix+7]+=i*(d*nbasis*ngGL*fpe);// bfd
	  iptrc[ix+8]+=i*(d*d*ngGL*fpe);     // JinvF
	  iptrc[ix+9]+=i*(d*ngGL*fpe);       // faceWeight
	  iptrc[ix+10]+=i*(nbasis*nbasis);             // mass 

	  iptrc[ix+11]+=i*(fpe*3*nfields*ngGL);  //faceFlux
	  iptrc[ix+12]+=i*fpe; 			   // cut2neigh & cut2face & cutoverset
          iptrc[ix+13]+=i*(maxseg*ngGL*d);	// OSFxn
          iptrc[ix+14]+=i*(maxseg*ngGL);	// OSFeID 
          iptrc[ix+15]+=i*(maxseg*ngGL*nbasis);	// OSFshpL, OSFshpR
          iptrc[ix+16]+=i*(maxseg*ngGL*3*nfields); // OSFflux
	}

	CUT_CELLS(x0,
		  x.data(),
		  xcut.data(), iptr.data(),
		  cut2e.data(),
		  d, etype, p, nelem, pc,
		  cut2face.data(),cut2neigh.data(),
		  elem2face, faces, iblank.data(),cutoverset.data(),imesh,ngGL);
	
	for(int i=0;i<necut;i++)
	  printf("cut elem %i: neigh = %i %i %i\n",
		 i,cut2neigh[3*i+0],cut2neigh[3*i+1],cut2neigh[3*i+2]);
      }
    };

    void mass_matrix(int imesh) {
      MASS_MATRIX(mass.data(),x.data(),
		  iptr.data(),d,etype,p,nelem,pc, imesh);
    }

    void cut_metrics(double x0, int imesh) {
      if (necut > 0) {
	COMPUTE_CUT_METRICS(x.data(),
			    JinvV.data(),
			    detJ.data(),
			    JinvF.data(),
			    iptr.data(),
			    d,etype,p,pc,pccut,
			    xcut.data(),
			    bvcut.data(),
			    bvdcut.data(),
			    JinvVcut.data(),
			    JinvFcut.data(),
			    detJcut.data(),
			    bfcutL.data(),
			    bfdcutL.data(),
			    bfcutR.data(),
			    bfdcutR.data(),
			    fwcut.data(),
			    iptrc.data(),
			    necut,
			    cut2e.data(),
			    cut2neigh.data(),imesh,
			    cutoverset.data());
	
	CUT_MASS_MATRIX(mass.data(),
			x.data(),
			JinvV.data(),
			iptr.data(),
			xcut.data(),
			detJcut.data(),
			iptrc.data(),
			d,etype,p,nelem,
			pc,pccut,necut,
			cut2e.data(),imesh);
      }
    };

    void initTimeStepping(int imesh) {
      /* compute some statistics of the mesh and report them */
        printf("\n#---------MESH %i------------\n",imesh);
	printf("#(nnodes, nelem, p)=(%d, %d, %d)\n",nnodes,nelem,p);
        printf("#nbasis=%d\n",nbasis);
	printf("#ndof=%d\n",nelem*nbasis);
	printf("#nfaces=%d\n",nfaces);
	printf("#totalArea=%f\n",total_area(detJ.data(),etype,p,d,nelem));
	printf("#ireg=%d\n",ireg);
	printf("#nsteps=%d\n",nsteps);
	printf("#Input parameters = ");
	output_params();
	printf("\n#--------------------------\n");
	printf("#%s\t%10s\t%14s\t%10s\n","step","l2","linf-loc","linf");
	ndof=nfields*nbasis*nelem;
    };

    void setupOverset(std::vector<int>& iptrB,
		      std::vector<int>& iptrcB,
		      std::vector<double>& xB,
		      std::vector<double>& xcutB,
		      std::vector<double>& JinvB,
		      std::vector<int>& cut2eB,
		      int necutB, int nelemB)
    {
    SETUP_OVERSET(cut2e.data(),
		  cut2eB.data(),
		  cutoverset.data(),
                  iptr.data(),
                  iptrB.data(),
                  iptrc.data(),
                  iptrcB.data(),
                  x.data(),
                  xB.data(),
		  xcut.data(),
		  xcutB.data(),
                  bfcutL.data(),
                  bfcutR.data(),
		  JinvV.data(),
		  JinvB.data(),
		  OSFnseg.data(),
		  OSFeID.data(),
		  OSFxn.data(),
		  OSFshpL.data(),
		  OSFshpR.data(),
                  d, etype, p, pc, pccut,
                  necut, necutB, nelemB);
    }
    
    void exchangeOverset(std::vector<double>& qB,
		         std::vector<int>& iptrB,
			 int imesh)
    {
      EXCHANGE_OVERSET(OSFflux.data(),
		       OSFshpL.data(),
		       OSFshpR.data(),
		       OSFnseg.data(),
		       OSFeID.data(),
		       cut2e.data(),
		       q.data(), 
		       qB.data(), 
		       cutoverset.data(),
		       iptrc.data(),
		       iptr.data(),
		       iptrB.data(),
		       necut,pccut,d,etype,p,pc,pde,imesh);
    } 

    void computeRHS(std::vector<double>& qsrc,int imesh) {
      COMPUTE_RHS(R.data(),
		  mass.data(),
		  bv.data(),
		  bvd.data(),
		  JinvV.data(),
		  detJ.data(),
		  bf.data(),
		  bfd.data(),
		  JinvF.data(),
		  faceWeight.data(),
		  fnorm.data(),
		  fflux.data(),
		  x.data(),
		  qsrc.data(),
		  elem2face,
		  iptr.data(),
		  iptf.data(),
		  faces,
		  pc,pf,pccut,pde,d,etype,p,nfaces,nelem,
		  bvcut.data(),
		  bvdcut.data(),
		  detJcut.data(),
		  bfcutL.data(),
		  bfcutR.data(),
		  fwcut.data(),
		  fcflux.data(),
		  OSFnseg.data(),
		  OSFeID.data(),
		  OSFxn.data(),
		  OSFshpL.data(),
		  OSFshpR.data(),
		  OSFflux.data(),
		  iptrc.data(),
		  necut,
		  cut2e.data(),
		  cut2face.data(),
		  cut2neigh.data(),
		  iblank.data(),
		  ireg,cutoverset.data(),imesh);
    }

    void update(std::vector<double>&qdest, std::vector<double>&qsrc, double dtfac)
    {
      UPDATE_DOFS(qdest.data(),dtfac,qsrc.data(),R.data(),ndof);
    };

    void rnorm(int &imax, double &rmax, double &rnorm, double dtfac)
    {
      rmax=rnorm=0;
      for(int i=0;i<ndof;i++) 
	{
	  if (fabs(R[i])> rmax) {
	    imax=i/nfields;
	    rmax=fabs(R[i]);
          }
	  rnorm+=(R[i]*R[i]);
	}
      rnorm=sqrt(rnorm/ndof);
      rmax*=(dtfac);
    };
    
    void output(int meshid, int istep) {
      OUTPUT_TECPLOT(meshid,istep,x.data(),q.data(),pc,iptr.data(),pde,d,etype,p,nelem);
    }
    
    int getNsteps() { return nsteps; };
    int getNsave()  { return nsave;  };
    int getNecut()  { return necut;  };
    double getDt()  { return dt;};
};
