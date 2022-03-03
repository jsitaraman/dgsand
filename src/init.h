void init_data(double *x,double *q, double *X, double *Q, int d, int nfields, int e, int p)
{
  int i,j,k,iflag;
  double **V; // vandermonde matrix
  double **Y; // rhs (coordinate locations and field)
  double u[d];
  int nbasis;

  /* number of basis for given element */
  nbasis=order2basis[e][p];
  /* form the Vandermonde matrix based on known location of 
     data from eloc */
  V=(double **)malloc(sizeof(double *)*nbasis);
  for(i=0;i<nbasis;i++) {
    V[i]=(double *)malloc(sizeof(double)*nbasis);
    for(j=0;j<nbasis;j++)
      {
	for(k=0;k<d;k++)
	  u[k]=eloc[e][p][i*d+k];
	V[i][j]=basis[e][j](u);
      }
  }
  /* form the RHS matrix from physical values of coordinates and
     field values */

  Y=(double **)malloc(sizeof(double *)*(nfields+d));
  for(i=0;i<nfields+d;i++) 
    Y[i]=(double *)malloc(sizeof(double)*nbasis);
  
  for(j=0;j<nbasis;j++)
    {
      for(i=0;i<d;i++)
	Y[i][j]=X[d*j+i];
        printf("X(%i) = %f\n",d*j+i,Y[i][j]);    
      for(i=d;i<nfields+d;i++)
	Y[i][j]=Q[(i-d)*nbasis+j];
    }
  /* solve V Y' = Y, Y=Y' on return from function */ 
  solvec(V,Y,&iflag,nbasis,nfields+d);
  
  for(i=0;i<d;i++)
    for(j=0;j<nbasis;j++)
      x[i*nbasis+j]=Y[i][j];
  for(i=d;i<nfields+d;i++)
    {
      for(j=0;j<nbasis;j++)
	{
	  q[(i-d)*nbasis+j]=Y[i][j];
	  //printf("%16.12f ",Y[i][j]);
	}
      //printf("\n");
    }
  
  for(i=0;i<nbasis;i++) free(V[i]);
  for(i=0;i<nfields+d;i++) free(Y[i]);
  free(V);
  free(Y);    
}  


void INIT_FIELDS(double *xcoord,int *e2n, 
		 double *Q, double *x, double *q, // data populated here
		 int *iptr, int pde, int e, int p,
		 int d, int nbasis, int itype, int nelem, int pc)
{
  int nvert=order2basis[e][p+(p==0)];
  int nfields=get_nfields[pde](d);  
  double X[d*nvert];
  int i,j,k,m,iQ,iq,ix,f;
  double r,s,qv[nfields];

  for(i=0;i<nelem;i++)
    {
      iq=iQ=iptr[pc*i];
      ix=iptr[pc*i+1];
      if (p > 0) {
	m=0;
	for(j=0;j<nvert;j++)
	  for(k=0;k<d;k++)
	  X[m++]=xcoord[d*e2n[nvert*i+j]+k];
        init_fields[pde](X,&(Q[iQ]),d,nbasis,itype);      
	init_data(&(x[ix]),&(q[iq]),X,&(Q[iQ]),d,nfields,e,p);
      }
      else {
	m=0;
	for(k=0;k<d;k++)
	  {
	    X[m]=0;
	    for(j=0;j<nvert;j++)
	      {
		x[ix+k*nvert+j]=xcoord[d*e2n[nvert*i+j]+k];
		X[m]+=x[ix+k*nvert+j];
	      }
	    X[m]/=nvert;
	    m++;
	  }
        init_fields[pde](X,&(q[iq]),d,nbasis,itype);      
      }
    }  
}
