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
  for(w=0;w<ngElem[e][p];w++)
    {
      for(j=0;j<d;j++)
	u[j]=gauss[e][g][(d+1)*w+j];
      //xx=yy=0;
      for(b=0;b<nbasis;b++)
	{
	  for(j=0;j<d;j++)
	    bd[b][j]=basis_d[e][b*d+j](u);
	  if (p > 0) bv[l++]=basis[e][b](u); // filled in as bv[nGL][nbasis]
	  //xx+=x[b]*bv[l-1];
	  //yy+=x[nbasis+b]*bv[l-1];
	}
      if (p==0) bv[l++]=1;
      //printf("(xx,yy)=%f %f\n",xx,yy);
      /*
	for(i=0;i<d;i++)
	for(b=0;b<nbasis;b++)
	printf("%f ", x[i*nbasis+b]);
	printf("\n");
      */
      for(i=0;i<d;i++)
	{
	  for(j=0;j<d;j++)
	    {
	      mat[i][j]=x[i*nbasis]*bd[0][j];
	      for(b=1;b<nbasis+(nbasis==1);b++)
		mat[i][j]+=x[i*nbasis+b]*bd[b][j];
	    }
	}
      
      if (d==2) invmat2x2(mat,jac,det);
      for(i=0;i<d;i++)
	for(j=0;j<d;j++) 
	  Jinv[ij++]=jac[i][j];
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
	  v=gaussgl[e][g][(d)*w];
	  wgt=gaussgl[e][g][(d)*w+1];	  

	  /* this is specific to 2-D */
	  /* remember fill order will be reversed for element
             sharing this edge */  
          f1=(f+1)%nfaces;
	  u[0]=(1-v)*eloc[e][p][d*f]+v*eloc[e][p][d*f1];
	  u[1]=(1-v)*eloc[e][p][d*f+1]+v*eloc[e][p][d*f1+1];
	  //xx=yy=0;
	  for(b=0;b<nbasis;b++)
	    {
	      for(j=0;j<d;j++)
		{
		  bd[b][j]=basis_d[e][b*d+j](u);
		  //bfd[ld++]=bd[b][j];      // filled in as bfd[nfaces][nGL][nbasis][d];
		}
	      if (p > 0) bf[l++]=basis[e][b](u);  // filled as bf[nfaces][nGL][nbasis]
	      //xx+=x[b]*bf[l-1];
	      //yy+=x[nbasis+b]*bf[l-1];
	    }
	  if (p==0) bf[l++]=1.0;
	  //printf("(r,s,xx,yy)=%f %f %f %f\n",r,s,xx,yy);
	  
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
                Ja[i]+=(mat[i][j]*face2elem[e][d*f+j]);
		//TODO add Jb[i] calculation for 3D elements here
	       }
            }
	    cross(&(faceWeight[2*m]),Ja,Jb,d);
            //printf("fw:(%f %f)\n",faceWeight[2*m],faceWeight[2*m+1]);
            if (d==2) invmat2x2(mat,jac,det);
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
			  int *iptr, int d, int e, int p, int nelem, int pc)
{
  int i,j,b;
  int ip,ix,ibv,ibvd,ibf,ij,idetj,ibfd,ijf,ifw;
  
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

      //Jacobian(&(x[ix]),&(bv[ibv]),&(bvd[ibvd]),&(JinvV[ij]), 
      //	       &(JinvBkD[ibk]), &(detJ[idetj]),d,e,p);

      Jacobian(x+ix, bv+ibv, bvd+ibvd, JinvV+ij,detJ+idetj,d,e,p);
      FaceWeights(x+ix,bf+ibf,bfd+ibfd,JinvF+ijf,faceWeight+ifw,d,e,p);

      /*      FaceWeights(&(x[ix]), &(bf[ibf]), &(bf[ibfd]), &(JinvF[ijf]), &(faceWeight[ifw]), 
	      d,e,p);*/

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


