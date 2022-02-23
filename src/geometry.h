void mass_matrix(double *M, double *x, int d, int e, int p)
{
  int i,j,ij,w,b,ii,jj;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];
  double mat[d][d],jac[d][d];
  double u[d],wgt,det;
  int g=p2g[e][p+1]; // gauss point type for this element type, use one
                     // order higher making sure mass matrix is exact

  if (p > 0 ) {
    for(i=0;i<nbasis;i++)
      for(j=0;j<nbasis;j++)
	{
	  ij=nbasis*i+j;
	  M[ij]=0;
	  for(w=0;w<ngElem[e][p+1];w++) // 
	    {
	      for(jj=0;jj<d;jj++)
		u[jj]=gauss[e][g][(d+1)*w+jj];
	      /* evaluate the jacobian, it's not stored at this many locations */
	      for(b=0;b<nbasis;b++)
		for(jj=0;jj<d;jj++)
		  bd[b][jj]=basis_d[e][b*d+jj](u);

	      for(ii=0;ii<d;ii++)
		{
		  for(jj=0;jj<d;jj++)
		    {
		      mat[ii][jj]=x[ii*nbasis]*bd[0][jj];
		      for(b=1;b<nbasis;b++)
			mat[ii][jj]+=x[ii*nbasis+b]*bd[b][jj];
		    }
		}
	      
	      if (d==2) invmat2x2(mat,jac,det);
	      wgt=gauss[e][g][(d+1)*w+2]*det;
	      /* could use bv here instead of reevaluating */
	      /* revaluate this if mesh is deforming */
	      M[ij]+=(wgt*basis[e][i](u)*basis[e][j](u));
	    }
	}
  }
  else {
    /* p=0 */
      for(b=0;b<nbasis;b++)
       for(jj=0;jj<d;jj++)
          bd[b][jj]=basis_d[e][b*d+jj](u);
      
      for(ii=0;ii<d;ii++)
	{
	  for(jj=0;jj<d;jj++)
	    {
	      mat[ii][jj]=x[ii*nbasis]*bd[0][jj];
	      for(b=1;b<nbasis;b++)
		mat[ii][jj]+=x[ii*nbasis+b]*bd[b][jj];
	    }
	}
      if (d==2) invmat2x2(mat,jac,det);
      M[0]=0.5*det;
  }
/*
printf("------------------\n");
for(i=0;i<order2basis[e][p];i++)
{
for(j=0;j<order2basis[e][p];j++)
  printf("%16.12f ",M[i*nbasis+j]);
printf("\n");
}
printf("------------------\n");
*/
}

void Jacobian(double *x,double *bv, double *bvd, double *Jinv, 
              double *detJ, int d, int e, int p)
{
  int b,w,i,j,l,ld,ij,m,n;
  double u[d];
  double mat[d][d],jac[d][d],det;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];  // basis derivative  
  int g=p2g[e][p];  // gauss quadrature data for this element type
  //double xx,yy;
  l=n=m=ld=ij=0;
  for(w=0;w<ngElem[e][p];w++) // loop over gauss pts
    {
      for(j=0;j<d;j++) // get rs coordinates of gauss pt
	u[j]=gauss[e][g][(d+1)*w+j];
      //xx=yy=0;
      for(b=0;b<nbasis;b++) // loop over bases
	{
	  for(j=0;j<d;j++) 
	    bd[b][j]=basis_d[e][b*d+j](u); // accumulate bases derivs at quad pt
	  if (p > 0) bv[l++]=basis[e][b](u); // filled in as bv[nGL][nbasis]
	}
      if (p==0) bv[l++]=1;
      //build jacobian dx/dr
      for(i=0;i<d;i++)
	{
	  for(j=0;j<d;j++)
	    {
	      mat[i][j]=x[i*nbasis]*bd[0][j];
	      for(b=1;b<nbasis+(nbasis==1);b++)
		mat[i][j]+=x[i*nbasis+b]*bd[b][j];
	    }
	}
      //invert jacobian (stored in jac) and get detJ
      if (d==2) invmat2x2(mat,jac,det);
      for(i=0;i<d;i++)
	for(j=0;j<d;j++) 
	  Jinv[ij++]=jac[i][j];

      // get basis derivs dN/dx = dN/dr*(dx/dr)^-1
      if (p > 0) {
      for(b=0;b<nbasis;b++)
	for(i=0;i<d;i++)
	  {
	     bvd[ld]=jac[0][i]*bd[b][0];
	      for(j=1;j<d;j++)
		bvd[ld]+=jac[j][i]*bd[b][j];
	     ld++;
	  }
      }
      else {
        for(i=0;i<d;i++) bvd[ld++]=0;
      }
      detJ[n++]=det;
    }
}

void cross(double *a,double *b,double *c,int d)
{
   a[0]=b[1]*c[2]-b[2]*c[1];
   a[1]=b[2]*c[0]-b[0]*c[2];
   if (d == 3) { 
   a[2]=b[0]*c[1]-b[1]*c[0];
   }
   return;
}
  
void FaceWeights(double *x, double *bf, double *bfd, double *Jinv, double *faceWeight, 
		 int d, int e, int p)
{
  int b,w,i,j,ij,l,ld,n,m,f,f1;
  double u[d],v,wgt;
  double mat[d][d],jac[d][d],det;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];  // basis derivative  
  int nfaces=facePerElem[e];
  int g=p2gf[e][p]; // gauss quadrature type for this element type
  double xx,yy;
  // ok, we are not going above 3D now
  double Ja[3];
  double Jb[3]={0,0,1};

  // for every face
  m=l=0;
  ld=ij=0;

  //printf("---------------\n");
  for(f=0;f<nfaces;f++)
    {
      // for every Gauss-point on this face
      for(w=0;w<ngGL[e][p];w++)
	{
	  v=gaussgl[e][g][(d)*w];  // gauss location
	  wgt=gaussgl[e][g][(d)*w+1];	 // gauss weight

	  /* this is specific to 2-D */
	  /* remember fill order will be reversed for element
             sharing this edge */  
          f1=(f+1)%nfaces;
	  // u is coordinate of gauss pt in rst
	  u[0]=(1-v)*eloc[e][p][d*f]  +v*eloc[e][p][d*f1]; //eloc is node location in rst; doing 1D interp?
	  u[1]=(1-v)*eloc[e][p][d*f+1]+v*eloc[e][p][d*f1+1];
	  //xx=yy=0;
	  
	  // get the bases and derivs at gauss legendre pts
	  for(b=0;b<nbasis;b++)
	    {
	      for(j=0;j<d;j++)
		{
		  bd[b][j]=basis_d[e][b*d+j](u); 
		  //bfd[ld++]=bd[b][j];      // filled in as bfd[nfaces][nGL][nbasis][d];
		}
	      // basis[][](u) computes basis function at u
	      if (p > 0) bf[l++]=basis[e][b](u);  // filled as bf[nfaces][nGL][nbasis]
	      //xx+=x[b]*bf[l-1];
	      //yy+=x[nbasis+b]*bf[l-1];
	    }
	  if (p==0) bf[l++]=1.0;
	  //printf("(r,s,xx,yy)=%f %f %f %f\n",r,s,xx,yy);

	  // build the jacobians of the edge
	  for(i=0;i<d;i++)
	    {
              faceWeight[d*m+i]=0;
	      Ja[i]=0;
	      for(j=0;j<d;j++)
		{
	          mat[i][j]=x[i*nbasis]*bd[0][j];
		  for(b=1;b<nbasis+(nbasis==1);b++)
		    mat[i][j]+=x[i*nbasis+b]*bd[b][j]; 
		} 
	      /* need a cross product here */
              for(j=0;j<d;j++)
	       {
                Ja[i]+=(mat[i][j]*face2elem[e][d*f+j]); // face2elem is...?
		//TODO add Jb[i] calculation for 3D elements here
	       }
            }
	    // do faceWeight = Ja x zhat
	    // This gives me the normal vector
	    cross(&(faceWeight[2*m]),Ja,Jb,d); 
            //printf("fw:(%f %f)\n",faceWeight[2*m],faceWeight[2*m+1]);

            // get [dx/dr]^-1
            if (d==2) invmat2x2(mat,jac,det);

            // compute [dr/dx][dN/dr]
	    if (p > 0) {
 	     for(b=0;b<nbasis;b++)
	      {
		for(i=0;i<d;i++)
		  {
		    bfd[ld]=jac[0][i]*bd[b][0];
		    for(j=1;j<d;j++)
		       bfd[ld]+=jac[j][i]*bd[b][j];
		    ld++;
		  }
	      }
	    }
	    else {
              for(i=0;i<d;i++) bfd[ld++]=0;
	    }
	    for(i=0;i<d;i++)
              for(j=0;j<d;j++)
	          Jinv[ij++]=jac[i][j];
	    m++;
	}
    }
}

void COMPUTE_GRID_METRICS(double *x, double *bv, double *bvd,double *JinvV, 
			  double *detJ,double *bf, double *bfd, double *JinvF, double *faceWeight,
			  int *iptr, int d, int e, int p, int nelem, int pc,
			  XXX)// need to add inputs for cut regions
{
  int i,j,b;
  int ip,ix,ibv,ibvd,ibf,ij,idetj,ibfd,ijf,ifw;
  
  // Metrics for the uncut cells
  for(i=0;i<nelem;i++)
    {
      ip=pc*i;
      ix   =iptr[ip+1];
      ibv  =iptr[ip+2];
      ibvd =iptr[ip+3];
      ij   =iptr[ip+4];
      idetj=iptr[ip+5];

      ibf  =iptr[ip+6];
      ibfd =iptr[ip+7];
      ijf  =iptr[ip+8];
      ifw  =iptr[ip+9];

      Jacobian(x+ix, bv+ibv, bvd+ibvd, JinvV+ij,detJ+idetj,d,e,p); // basis on vol
      FaceWeights(x+ix,bf+ibf,bfd+ibfd,JinvF+ijf,faceWeight+ifw,d,e,p); // basis on face
    }

  // Metrics for the cut regions
  for(i=0;i<necut;i++)
  {
    ip = eid[i]*pc; // Get elem ID of original elem somehow XXX

    // Quantities for the original element
    ix   =iptr[ip+1]; 
    ibvd =iptr[ip+3];
    ij   =iptr[ip+4];
    idetj=iptr[ip+5];

    ibf  =iptr[ip+6];
    ibfd =iptr[ip+7];
    ijf  =iptr[ip+8];
    ifw  =iptr[ip+9];

    // Quantities for the cut element
    cip = pc*i;
    ecut = xxx; // array of cut cell element type?

    cix   =iptrc[cip+1]; // x coord
    cibvd =iptrc[cip+3]; // shp func deriv
    cij   =iptrc[cip+4]; // Jinv
    cidetj=iptrc[cip+5]; // detJ

    cibf  =iptrc[cip+6]; // face shp function
    cibfd =iptrc[cip+7]; // face shp deriv 
    cijf  =iptrc[cip+8]; // JinvF
    cifw  =iptrc[cip+9]; // faceWeight

    BasesVCut(XXX); // get bases at cut vol quad pts
    FaceWeights(x+cix,bf+cibf,bfd+cibfd,JinvF+cijf,faceWeight+cifw,d,ecut,pcut); // basis on face
  }

}

void MASS_MATRIX(double *mass,double *x, int *iptr, int d, int e, int p, int nelem, int pc)
{
  int i;
  int ix,im;
  for(i=0;i<nelem;i++)
    {
      ix=iptr[pc*i+1];
      im=iptr[pc*i+10];
      mass_matrix(&(mass[im]),&(x[ix]),d,e,p);
    }
}

void BasesVCut(XXX)
// Routine that builds up the bases and derivs of the interior quad pts
{
  //Loop over vol quad pts  
  double rstcut[d]; 
  for(w=0; w<ngElem[ecut][pcut]; w++){

    // get local coord of quad pt
    for(j=0; j<d;j++) rstcut[j] = gauss[ecut][g][(d+1)*w+j];

    // get bases and derivs at quad pt
    // XXX double check that the inputs are correct
    CutCellInterp(x,d,e,Jinv,rstcut,xcut,ecut,pcut,bvcut[xxx],bvdcut[xxx],invJcut,detJcut);
  }

}

void CutCellInterp(double *x, int d, int e, double* Jinv, 
   		   double *rstcut, double *xcut, int ecut, int pcut, 
		   double* bvcut, double*bvdcut, double *invJcut, double *detJcut)
//
// Routine that outputs bases and basis derivs at a specified rst 
// location on the cut cell. Used for both face and volume quad pts
//
// Inputs: 
//   x = global coords of orig tri
//   d = number of spat dims
//   e = original element type
//   Jinv = orig element inv Jac
//   rstcut = cut cell rst value of desired pt
//   xcut = global coords of cut region
//   ecut = cut region element type
//   pcut = cut region order of integration (same as orig p)
//
// Outputs
//   bcut = bases values at cut cell gauss pts
//   bdcut = bases derivs at cut cell gauss pts
//   invJcut = inv Jac of cut cell
//   detJcut = det J of cut cell
{
 
  // notation: 
  //   xyz = global coord sys
  //   rst = local coord sys of original cell
  //   ijk = local coord sys of cut cell

  int b,w,i,j,l,ld,ij,m,n;
  int nbasis; // number of basis functions
  double u[d],jcut[d*d],mat[d][d],ijac[d][d],det;
  double xycut[d*d];
  int g=p2g[ecut][pcut];
  
  l=n=m=ld=ij=0;

  //get xyz coord of nodes on the two elements
  //corresponding to rst = 0,0,1
  double x0 = x[XXX]; 
  double x0cut = xcut[XXX];
 
  nbasis=order2basis[ecut][1]; // XXX double check last argument, should be p=1
  for(b=0;b<nbasis;b++){
    for(j=0;j<d;j++) bd[b][j]=basis_d[ecut][b*d+j](rstcut); // accumulate bases derivs at quad pt
  } // loop over bases
        
  // Compute p=1 jacobian for subcell element, dx/di
  // only p=1 required for this regardless of element p 
  for(i=0;i<d;i++){
    for(j=0;j<d;j++){
      mat[i][j]=xcut[i*nbasis]*bd[0][j];
      for(b=1;b<nbasis+(nbasis==1);b++)
        mat[i][j]+=x[i*nbasis+b]*bd[b][j];
    } // loop over d
  } // loop over d

  // Get and store inverse and detJ for cut cell
  if(d==2) invmat2x2(mat,ijac,det)
  detJcut[n++]=det; 

  // reshape
  for(i=0;i<d;i++)
    for(j=0;j<d;j++){
      jcut[(i-1)*d+j] = mat[i][j]; 
      invJcut[ij++]=ijac[i][j];
    }

  // Compute the global x y location of this cut cell quad pt
  // [xycut] = [jac_cut]*[ij] + x0_cut
  axb(jcut,u,xycut,d);
  for(j=0;j<d;j++) xycut[d] = xycut[d] + x0cut[d] - x0[d];
    
  // Compute rst coord of cut cell quad
  // [rs] = [Jinv_orig][xycut - x0_orig]
  axb(Jinv,xycut,u,d);

  // Store orig cell bases value and deriv at cut cell quad pt
  // XXX double check this. is Jinv in correct order?
  if(p>0){
    for(b=0;b<nbasis;b++){
      for(i=0;i<d;i++){
        bvcut[l++]=basis[e][b](u); // filled in as b[nGL][nbasis]
        for(j=0;j<d;j++) bvdcut[ld]=Jinv[xxx]*basis_d[e][b*d+j](u);
        ld++; 
      }
    }
  }
}

void axb(double* a, double* x, double* b, int d)
{

  int ind; 
  for(int i=0;i<d;i++){
    b[i]=0; 
    for(int j=0;j<d;j++){
      ind = i*d+j; 
      b[i] += a[ind]*x[j];
    }
  }

}
