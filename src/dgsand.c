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
#include "inputParser.h"
#include "solvec.h"
#include "basislib_u.h"
#include "quadrature.h"
#include "memutil.h"
#include "pde.h"
#include "find_faces.h"
#include "geometry.h"
#include "init.h"
#include "residual.h"
#include "read_grid.h"
#include "output.h"

void main(void)
{
  double *q,*qstar,*Q,*R;                           // field
  double *x,*bv,*bvd,*JinvV,*detJ;  // volumetric quantities
  double *xcut,*bvcut,*bvdcut,*JinvVcut,*detJcut;  // volumetric quantities
  double *bf,*bfd,*JinvF,*faceWeight;        // face quantities
  double *bfcut,*bfdcut,*JinvFcut,*fwcut;        // face quantities
  double *mass;                              // mass matrix
  double *fnorm,*fflux;                      // face normals and face flux
  double *fcnorm,*fcflux;                      // face normals and face flux
  double xsum;

  /* inputs */
  int pde;                          // pde type, only Euler eqns are implemented now
  int d=2;                          // dimensions (only implemented for 2D now)
  int nfields;                      // number of fields for this pde
  int p;                            // p order
  int etype=0;                      // element type (triangle = 0, quad = 1)
  
  int nbasis;                // number of bases for solution
  int nbasisx;                      // number of bases for grid
  int itype;                        // initialization type
  
  /* overset inputs */
  int nmesh = 2; 

  /* input and post processed from a grid file */
  /* Hacking this and assuming that both grids are identical */
  int nnodes;                       // number of nodes in grid
  int nelem;                        // number of elements in each mesh
  int nbnodes;                      // number of primal nodes on physical boundaries
  int *ibc;                         // boundary condition node indices and their types
  int nfaces;                       // number of faces 
  double *xcoord;                   // coordinates of provided grid (contains both grids)
  int *elem2node;             // element to node connectivity
  int *faces,*elem2face;     // face to cell connectivity and element to face information
  int *iptr,*iptf;                  // pointer into data arrays
  int *iptrc,*iptrcf;                  // pointer into data arrays


  /* local variables */
  int a,i,j,m,b,ix, pc,pf,fpe,imax,ndof,n,nsteps,nsave;
  double wgt,totalArea,rnorm,rmax,dt;
  /* rk3 coefficients */
  double rk[4]={0.25,8./15,5./12,3./4};

  /* parse inputs */
  parseInputs("input.dgsand",&pde,&itype,&nsteps,&dt,&nsave);
  nfields=get_nfields[pde](d);

  /* read a 2D grid */
  readgrid2D(&xcoord,&elem2node,&ibc,&p,&nnodes,&nelem,&nbnodes);
  nbasis=order2basis[etype][p];         // basis for solution
  nbasisx=order2basis[etype][p+(p==0)]; // basis for grid

  /* find element to face connectivity */
  find_faces(elem2node,&elem2face,&faces,&nfaces,ibc,nelem,nbnodes,3,nbasisx);

  /* allocate memory */
  /* field parameters per element */  
  q=dgsand_alloc(double,(nfields*nbasis*nelem));     // modal coefficients
  qstar=dgsand_alloc(double,(nfields*nbasis*nelem)); // modal coefficients
  Q=dgsand_alloc(double,(nfields*nbasis*nelem));     // values at physical locations
  R=dgsand_alloc(double,(nfields*nbasis*nelem));     //solution residual 

  /* geometrical parameters per volume QP of each element */
  /* TODO: some of these such as bv and JinvV can be optimized or omitted */
  x       =dgsand_alloc(double,(d*(nbasisx))*nelem);              // coord modal coefficients (p0 mod)
  bv      =dgsand_alloc(double,(nbasis*ngElem[etype][p]));        // basis value at volume QP
  bvd     =dgsand_alloc(double,(d*nbasis*ngElem[etype][p]*nelem));// basis derivative value at volume QP
  JinvV   =dgsand_alloc(double,(d*d*ngElem[etype][p]*nelem));     // J^{-1} at volume QP
  detJ    =dgsand_calloc(double,(ngElem[etype][p]*nelem));        // |J| at volume QP

  /* geometrical parameters per face QP of each element */
  /* TODO : some these such as bf and JinvF can optimized/omitted */
  fpe = facePerElem[etype];
  bf        =dgsand_alloc(double,(nbasis*ngGL[etype][p]*fpe));  // basis value at face QP
  bfd       =dgsand_alloc(double,(d*nbasis*ngGL[etype][p]*fpe*nelem));// basis der. value at face QP
  JinvF     =dgsand_alloc(double,(d*d*ngGL[etype][p]*fpe*nelem));     // J^{-1} at face QP
  faceWeight=dgsand_alloc(double,(d*ngGL[etype][p]*fpe*nelem));       // faceNormals at face QP
  mass      =dgsand_alloc(double,(nbasis*nbasis*nelem));              // mass matrix

  fnorm     =dgsand_alloc(double,(d*ngGL[etype][p]*nfaces));          // face normals
  fflux     =dgsand_alloc(double,(3*nfields*ngGL[etype][p]*nfaces));  // face fields and flux        

  /* pointer array into each data array above */
  pc=11; // number of unique sizes with elements
  pf=2;  // number of unique sizes associated with faces
  iptr=dgsand_calloc(int,(pc*nelem));
  iptf=dgsand_calloc(int,(pf*nfaces));
 
  /* set the pointers, TODO: this has to change when there is a variety of elements */
  for(i=0;i<nelem;i++){    
      ix=a*pc*nelem+pc*i;
      iptr[ix]+=i*(nfields*nbasis);               // q, Q, R
      iptr[ix+1]+=i*(d*(nbasisx));                // x
      iptr[ix+2]+=0;                              // bv (this is same per element type)
      iptr[ix+3]+=i*(d*nbasis*ngElem[etype][p]);  // bvd
      iptr[ix+4]+=i*(d*d*ngElem[etype][p]);       // JinvV
      iptr[ix+5]+=i*(ngElem[etype][p]);           // detJ

      iptr[ix+6]+=0;                              // bf (this is same per element type)
      iptr[ix+7]+=i*(d*nbasis*ngGL[etype][p]*fpe);// bfd
      iptr[ix+8]+=i*(d*d*ngGL[etype][p]*fpe);     // JinvF
      iptr[ix+9]+=i*(d*ngGL[etype][p]*fpe);       // faceWeight
      iptr[ix+10]+=i*(nbasis*nbasis);             // mass 
  }
  for(i=0;i<nfaces;i++){
      ix=a*pf*nfaces+pf*i;
      iptf[ix]+=(i*d*ngGL[etype][p]);            //faceNormal
      iptf[ix+1]+=(i*3*nfields*ngGL[etype][p]);  //faceFlux
  }

  /* Cut region parameters */ 
  // XXX need to decide how to handle etype here.
  // will we need an array to tell us how many nodes are in each element
  // or will we keep all cut cells as triangles?
  //
  // Assume for the hacked problem that the overlap only happens in the boundary elements
  //
  int necut = 0; 
  int ncfaces = 0; 

  xcut       =dgsand_alloc(double,(d*(nbasisx))*necut);              // coord modal coefficients (p0 mod)
  bvcut      =dgsand_alloc(double,(necut*nbasis*ngElem[etype][p]));        // basis value at volume QP
  bvdcut     =dgsand_alloc(double,(necut*d*nbasis*ngElem[etype][p]));// basis derivative value at volume QP
  JinvVcut   =dgsand_alloc(double,(necut*d*d*ngElem[etype][p]));     // J^{-1} at volume QP
  detJcut    =dgsand_calloc(double,(necut*ngElem[etype][p]));        // |J| at volume QP
  
  // Need to think more on these
  // Are they the right size? XXX
  fpe = facePerElem[etype];
  ncfaces = fpe*necut; // XXX is this right? 
  bfcut     =dgsand_alloc(double,(nbasis*ngGL[etype][p]*fpe*necut));  // basis value at face QP
  bfdcut    =dgsand_alloc(double,(d*nbasis*ngGL[etype][p]*fpe*necut));// basis der. value at face QP
  JinvFcut  =dgsand_alloc(double,(d*d*ngGL[etype][p]*fpe*necut));     // J^{-1} at face QP
  fwcut     =dgsand_alloc(double,(d*ngGL[etype][p]*fpe*necut));       // faceNormals at face QP

  fcnorm     =dgsand_alloc(double,(d*ngGL[etype][p]*ncfaces));          // face normals
  fcflux     =dgsand_alloc(double,(3*nfields*ngGL[etype][p]*ncfaces));  // face fields and flux        

  /* pointer array into each data array above */
  iptrc=dgsand_calloc(int,(pc*necut));
  iptrcf=dgsand_calloc(int,(pf*ncfaces));
  
  // XXX how do I handle fpe and element type being different for certain cut regions?
  // Assume for now all cuts are triangles?
  for(i=0;i<necut;i++){
      ix=a*pc*necut + pc*i;
      iptrc[ix]+=i*(nfields*nbasis);               // q, Q, R
      iptrc[ix+1]+=i*(d*(nbasisx));                // x
      iptrc[ix+2]+=i*(nbasis*ngElem[etype][p]);   // bv (this is NOT same per element type)
      iptrc[ix+3]+=i*(d*nbasis*ngElem[etype][p]);  // bvd
      iptrc[ix+4]+=i*(d*d*ngElem[etype][p]);       // JinvV
      iptrc[ix+5]+=i*(ngElem[etype][p]);           // detJ

      iptrc[ix+6]+=i*(nbasis*ngGL[etype][p]*fpe); // bf (this is NOT same per element type)
      iptrc[ix+7]+=i*(d*nbasis*ngGL[etype][p]*fpe);// bfd
      iptrc[ix+8]+=i*(d*d*ngGL[etype][p]*fpe);     // JinvF
      iptrc[ix+9]+=i*(d*ngGL[etype][p]*fpe);       // faceWeight
      iptrc[ix+10]+=i*(nbasis*nbasis);             // mass 
  }
  for(i=0;i<ncfaces;i++){
      ix=a*pf*ncfaces + pf*i;
      iptrcf[ix]+=(i*d*ngGL[etype][p]);            //faceNormal
      iptrcf[ix+1]+=(i*3*nfields*ngGL[etype][p]);  //faceFlux
  }
 
  /* initialize fields on all the elements */
  INIT_FIELDS(xcoord,elem2node,Q,x,q,iptr,pde,etype,p,d,nbasis,itype,nelem,pc);

  /* compute grid metrics */
  COMPUTE_GRID_METRICS(x,bv,bvd,JinvV,detJ,
  		       bf,bfd,JinvF,faceWeight,iptr,d,etype,p,nelem,pc);

  /* Do the blanking here based on detJ */

  /* compute the mass matrix for each element */
  MASS_MATRIX(mass,x,iptr,d,etype,p,nelem,pc);


  /* Handle cut cells */
/*
  if(necut>0){
    COMPUTE_CUT_METRICS(x,JinvV,detJ,
                        JinvF,iptr,d,etype,p,pc,
                        xcut,bvcut,bvdcut,JinvVcut,detJcut,
                        bfcut,bfdcut,JinvFcut,fwcut,
                        iptrc,necut);

    CUT_MASS_MATRIX(mass,x,iptr,d,etype,p,nelem,pc,necut);
  }
*/


  /* compute some statistics of the mesh and report them */
  totalArea=0.0;
  for(i=0;i<nelem;i++) 
    for(j=0;j<ngElem[etype][p];j++)
      {
	wgt=0.5*gauss[etype][p2g[etype][p]][(d+1)*j+2];
	totalArea+=(wgt*detJ[m]);
	m++;
      }
  printf("#---------dgsand------------\n");
  printf("#(nnodes, nelem, p)=(%d, %d, %d)\n",nnodes,nelem,p);
  printf("#ndof=%d\n",nelem*nbasis);
  printf("#nfaces=%d\n",nfaces);
  printf("#totalArea=%f\n",totalArea);
  printf("#Input parameters = ");
  for(i=0;i<6;i++) printf("%f ",param[i]);
  printf("\n#--------------------------\n");
  printf("#%s\t%10s\t%14s\t%10s\n","step","l2","linf-loc","linf"); 
  OUTPUT_TECPLOT(0,x,q,pc,iptr,pde,d,etype,p,nelem);

  ndof=nfields*nbasis*nelem;

  /* basic 3rd order RK time stepping */
  
  for(n=1;n<=nsteps;n++)
    {
      // XXX need to add inputs for cut cells
      COMPUTE_RHS(R,mass,bv,bvd,JinvV,detJ,
		  bf,bfd,JinvF,faceWeight,fnorm,fflux,
		  x,q,elem2face,iptr,iptf,faces,
		  pc,pf,pde,d,etype,p,nfaces,nelem,
                  bvcut,bvdcut,JinvVcut,detJcut,
                  bfcut,bfdcut,JinvFcut,fwcut,fcnorm,fcflux,
                  xcut,iptrc,iptrcf,necut);
      
      UPDATE_DOFS(qstar,rk[1]*dt,q,R,ndof);
      UPDATE_DOFS(q,rk[0]*dt,q,R,ndof);
      
      COMPUTE_RHS(R,mass,bv,bvd,JinvV,detJ,
		  bf,bfd,JinvF,faceWeight,fnorm,fflux,
		  x,qstar,elem2face,iptr,iptf,faces,
		  pc,pf,pde,d,etype,p,nfaces,nelem,
                  bvcut,bvdcut,JinvVcut,detJcut,
                  bfcut,bfdcut,JinvFcut,fwcut,fcnorm,fcflux,
                  xcut,iptrc,iptrcf,necut);
      
      UPDATE_DOFS(qstar,rk[2]*dt,q,R,ndof);

      COMPUTE_RHS(R,mass,bv,bvd,JinvV,detJ,
		  bf,bfd,JinvF,faceWeight,fnorm,fflux,
		  x,qstar,elem2face,iptr,iptf,faces,
		  pc,pf,pde,d,etype,p,nfaces,nelem,
                  bvcut,bvdcut,JinvVcut,detJcut,
                  bfcut,bfdcut,JinvFcut,fwcut,fcnorm,fcflux,
                  xcut,iptrc,iptrcf,necut);


      UPDATE_DOFS(q,rk[3]*dt,q,R,ndof);
      
      rmax=rnorm=0.0;
      for(i=0;i<ndof;i++) 
      {
         if (fabs(R[i])> rmax) {
          imax=i/nfields;
          rmax=fabs(R[i]);
          }
	 rnorm+=(R[i]*R[i]);
      }
      rnorm=sqrt(rnorm/ndof);
      rmax*=(rk[3]*dt);
      printf("%d\t%18.16f\t%d\t%18.16f\n",n,rnorm,imax,rmax);
      if (n%nsave==0) OUTPUT_TECPLOT(n,x,q,pc,iptr,pde,d,etype,p,nelem);
    }
}
