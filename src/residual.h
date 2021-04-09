
void volIntegral(double *residual, double *bv, double *bvd, double *q, double *detJ,
		 int pde, int d, int e, int p)
{
  int b,w,i,j,l,ld,m,f;
  int nbasis=order2basis[e][p];
  double wgt;
  double bvv[nbasis];
  double bvvd[nbasis][d];
  int nfields=get_nfields[pde](d);
  double flux[d][nfields];
  double qv[nfields],qvd[nfields][d];
  int g=p2g[e][p];
  
  l=ld=0;
  // for all gauss-quadrature points 
  for(w=0;w<ngElem[e][p];w++)
    {
      wgt=gauss[e][g][(d+1)*w+2]*detJ[w];	  
      /* collect basis values for this quadrature point */
      for(b=0;b<nbasis;b++) {
  	bvv[b]=bv[l++];
	for(j=0;j<d;j++)
	  bvvd[b][j]=bvd[ld++];
      }
 
      /* project the q field to the gauss point */
      for(f=0;f<nfields;f++)
	{
	  qv[f]=bvv[0]*q[f*nbasis];
	  for(b=1;b<nbasis;b++)	    
	    qv[f]+=(bvv[b]*q[f*nbasis+b]);
	  for(j=0;j<d;j++)
	    {
	      qvd[f][j]=bvvd[0][j]*q[f*nbasis];
	      for(b=1;b<nbasis;b++)
		qvd[f][j]+=(bvvd[b][j]*q[f*nbasis+b]);
	    }
	}
      /* compute the flux function in each dimension */
      for(j=0;j<d;j++)
	flux_function[pde](flux[j],qv,qvd,j);
      m=0;
      for(f=0;f<nfields;f++)
	for(b=0;b<nbasis;b++)
	  {
	    if (w==0) residual[m]=0;
	    for(j=0;j<d;j++)
	      residual[m]+=wgt*(bvvd[b][j]*flux[j][f]); /* \grad b . F */
	    m++;
	  }
    }
//  printf("resv :");
//  for(b=0;b<nbasis;b++) printf(" %f ",residual[b]);
//  printf("\n");
}


void faceIntegral(double *residual, double *fflux, double *bf, double *bfd, int *elem2face,
		  int *iptrf,double *q, int pf, int pde, int d, int e, int p, int ielem)
{
  int b,w,i,j,l,ld,m,f,fid,fst,fsgn;
  int nbasis=order2basis[e][p];
  double wgt,v;
  double bvv[nbasis];
  double bvvd[nbasis][d];
  int nfields=get_nfields[pde](d);
  double flux[d][nfields];
  double qv[nfields],qvd[nfields][d];
  int g=p2gf[e][p];
  int nfp=facePerElem[e];
  double resf[nfields*nbasis];

  l=ld=m=0;
  for(f=0;f<nfields;f++)
    for(b=0;b<nbasis;b++)
        resf[m++]=0; 

  for(i=0;i<nfp;i++)
    {
      fsgn=elem2face[i]/abs(elem2face[i]);
      fid=abs(elem2face[i])-1;
      //fst=iptrf[pf*fid+1]+2*nfields;
      fst=iptrf[pf*(fid+(1-fsgn)/2)+1]-(1-fsgn)*nfields/2+(1+fsgn)*nfields;
      //if (ielem==7157) printf("%d %d ",fsgn,fid);
      for(w=0;w<ngGL[e][p];w++)
	{
	  v=gaussgl[e][g][(d)*w];
	  wgt=gaussgl[e][g][(d)*w+1];
	  for(b=0;b<nbasis;b++) {
	    bvv[b]=bf[l++];
	    for(j=0;j<d;j++)
	      bvvd[b][j]=bfd[ld++];
	  }

	  m=0;
	  for(f=0;f<nfields;f++)
	    {
              //if (f==0) printf("%f %d ",fflux[fst+f],fsgn);
	      for(b=0;b<nbasis;b++)
		{
                  //if (f==0) printf("%f %f ",wgt,bvv[b]);
                  resf[m]-=(wgt*fsgn*fflux[fst+f]*bvv[b]);
		  //if (ielem==7157) printf("%f ",fflux[fst+f]);
		  residual[m]-=(wgt*fsgn*fflux[fst+f]*bvv[b]);
		  m++;
		}
	      //if (f==0) printf("\n");
	    }
	  fst+=(3*fsgn*nfields);
	}
      //if (ielem==7157) printf("\n");
    }
  //if (ielem==7157) {
  //printf("resf :");
  //for(f=0;f<nfields;f++)
  // for(b=0;b<nbasis;b++) printf(" %f ",residual[f*nbasis+b]);
  //printf("\n");
  //}
}

void setFaceQuantities(double *fnorm,double *fflux,int *elem2face, int *iptrf,
		       double *faceWeight, double *bf,double *bfd, double *q, 
		       int nfields, int pf, int pde, int d , int e, int p)
{
  int b,w,i,j,k,l,f,n,fid,floc,fst,fsgn,nst;
  int nbasis=order2basis[e][p];
  double bvv[nbasis],qv[nfields];
  int nfp=facePerElem[e];
  l=k=0;
  for(i=0;i<nfp;i++)
    {
      fsgn=elem2face[i]/abs(elem2face[i]);
      fid=abs(elem2face[i])-1;
      // pick out the right location for inserting fields for this face
      // the cell with negative sign for the edge fills in backward order
      nst=iptrf[pf*fid];
      fst=iptrf[pf*(fid+(1-fsgn)/2)+1]-(1-fsgn)*nfields;
      n=0;
      for(w=0;w<ngGL[e][p];w++)
	{
	  for(b=0;b<nbasis;b++)
	    bvv[b]=bf[l++];
	  for(f=0;f<nfields;f++)
	    {
	      floc=fst+f;
	      fflux[floc]=0; //bvv[0]*q[f*nbasis];
	      for(b=0;b<nbasis;b++)	    
		fflux[floc]+=(bvv[b]*q[f*nbasis+b]);
	    }	  
	  for(j=0;j<d*fsgn;j++)
	    fnorm[nst+(n++)]=faceWeight[k+j];
          k+=d;
	  fst+=(3*fsgn*nfields);
	}
    }
}

//void setFaceQuantities(double *fnorm,double *fflux,int *elem2face, int *iptrf,
//		       double *faceWeight, double *bf,double *bfd, double *q, 
//		       int nfields, int pde, int d , int e, int p)
void FILL_FACES(double *fnorm, double *fflux, int *elem2face,int *iptr, int *iptrf,
		double *faceWeight, double *bf, double *bfd, double *q, 
		int pc, int pf, int pde, int d, int e, int p, int nelem, int nfaces)
{
  int i,ix;
  int ifw,ibf,ibfd,iq;
  int nfields=get_nfields[pde](d);
  int nfp=facePerElem[e];
  double *dummy;
  //dummy=(double *)malloc(sizeof(double) *1000);
  for(i=0;i<nelem;i++)
    {
      ix=pc*i;
      iq=iptr[ix]; 
      ibf=iptr[ix+6];
      ibfd=iptr[ix+7];
      ifw=iptr[ix+9];
      
      setFaceQuantities(fnorm,fflux,elem2face+nfp*i, iptrf, 
			faceWeight+ifw, bf+ibf, bfd+ibfd, q+iq,
			nfields, pf, pde, d, e, p);
    }
}

void FILL_BC(double *fnorm,double *fflux, int *faces,
	     int pde,int d, int e, int p, int nfaces)
{
  int i,j,ifp;
  int nfields=get_nfields[pde](d);
  for(i=0;i<nfaces;i++)
    {
      /* TODO                                        */
      /* there should be a way to avoid this if loop */
      /* collect boundary faces in to a separate list in preproc */
      if (faces[6*i+4]==-1) {
	for(j=0;j<ngGL[e][p];j++)
	  {
	    ifp=(i*ngGL[e][p]+j)*3*nfields+nfields;
	    far_field[pde](fflux+ifp);
	  }
      }
    }
}
	     
	     
void COMPUTE_FACE_FLUXES(double *fnorm, double *fflux,
			 int pde, int d, int e, int p, int nfaces, int *faces)
{
  int nfields=get_nfields[pde](d);
  double normal[d],xnorm[d];
  int i,j,ifl,ifr,iflux,f;

  for(i=0;i<nfaces*ngGL[e][p];i++)
    {
      f=i/ngGL[e][p];
      for(j=0;j<d;j++)
	xnorm[j]=fnorm[d*i+j];

      ifl=i*3*nfields;
      ifr=ifl+nfields;
      iflux=ifr+nfields;

      gradient_indep_flux[pde](fflux+ifl,fflux+ifr,fflux+iflux,xnorm,0.0);
    }
}

void invertMass(double *mass, double *R, int pde, int d , int e, int p,int ielem)
{
  int i,j,f,b;
  int nbasis=order2basis[e][p];
  int iflag;
  int nfields=get_nfields[pde](d);
  /*
  if (ielem==7157) { 
  for(f=0;f<nfields;f++)
    {
      for(b=0;b<nbasis;b++)
	printf("%f ",R[f*nbasis+b]);
      printf("\n");
    }
  printf("mass=%f\n",mass[0]);
  printf("--------------------------\n");
  }   
  */
  solvec_copy_reshape(mass,R,&iflag,nbasis,nfields);
  /*
  //printf("iflag=%d\n",iflag);
  if (ielem==7157) {  
  for(f=0;f<nfields;f++)
    {
      for(b=0;b<nbasis;b++)
	  printf("%f ",R[f*nbasis+b]);
      printf("\n");
    }
  }
  */
  //for(i=0;i<nbasis;i++) free(mtmp[i]);
  //free(mtmp);
}


void checkGradients(double *x, double *q, double *bv, double *bvd, double *bf,
		    double *bfd, int pde, int d, int e, int p)
{

  int b,w,i,j,l,ld,m,f;
  int nbasis=order2basis[e][p];
  double wgt;
  double bvv[nbasis];
  double bvvd[nbasis][d];
  int nfields=get_nfields[pde](d);
  double flux[d][nfields];
  double xv[d],qv[nfields],qvd[nfields][d];
  int g=p2g[e][p];
  int nfp=facePerElem[e];
  
  // for all gauss-quadrature points
  l=ld=0;
  for(w=0;w<ngElem[e][p];w++)
    {
      wgt=gauss[e][g][(d+1)*w+2];	  
      /* collect basis values for this quadrature point */
      for(b=0;b<nbasis;b++) {
	bvv[b]=bv[l++];
	for(j=0;j<d;j++)
	  bvvd[b][j]=bvd[ld++];
      }
      /* project the q field to the gauss point */
      for(f=0;f<d;f++)
	{
	  xv[f]=bvv[0]*x[f*nbasis];
	  for(b=1;b<nbasis;b++)
	    xv[f]+=(bvv[b]*x[f*nbasis+b]);
	  //printf("%12.8f ",xv[f]);
	}
      for(f=0;f<nfields;f++)
	{
	  qv[f]=bvv[0]*q[f*nbasis];
	  for(b=1;b<nbasis;b++)	    
	    qv[f]+=(bvv[b]*q[f*nbasis+b]);
	  //printf("%12.8f ",qv[f]);
	  for(j=0;j<d;j++)
	    {
	      qvd[f][j]=bvvd[0][j]*q[f*nbasis];
	      for(b=1;b<nbasis;b++)
		qvd[f][j]+=(bvvd[b][j]*q[f*nbasis+b]);
	      //printf("%12.8f ",qvd[f][j]);
	    }
	}
      //printf("\n");
    }

  l=ld=0;
  for(i=0;i<nfp;i++)
    {
      // for every Gauss-point on this face
      for(w=0;w<ngGL[e][p];w++)
	{
	  for(b=0;b<nbasis;b++) {
	    bvv[b]=bf[l++];
	    for(j=0;j<d;j++)
	      bvvd[b][j]=bfd[ld++];
	  }	  
	  /* project the q field to the gauss point */
	  for(f=0;f<d;f++)
	    {
	      xv[f]=bvv[0]*x[f*nbasis];
	      for(b=1;b<nbasis;b++)
		xv[f]+=(bvv[b]*x[f*nbasis+b]);
	      printf("%12.8f ",xv[f]);
	    }
	  for(f=0;f<nfields;f++)
	    {
	      qv[f]=bvv[0]*q[f*nbasis];
	      for(b=1;b<nbasis;b++)	    
		qv[f]+=(bvv[b]*q[f*nbasis+b]);
	      printf("%12.8f ",qv[f]);
	      for(j=0;j<d;j++)
		{
		  qvd[f][j]=0; //bvvd[0][j]*q[f*nbasis];
		  for(b=0;b<nbasis;b++)
		    {
		      qvd[f][j]+=(bvvd[b][j]*q[f*nbasis+b]);
		    }
		   printf("%16.12f ",qvd[f][j]);
		}
	    }
	  printf("\n");
	}
    }

}
void CHECK_GRADIENTS(double *x, double *q, double *bv, double *bvd, double *bf,
		     double *bfd, int *iptr, int pc,
		     int pde, int d, int e, int p, int nelem)
{
  int ix,iq,ibv,ibvd,ibf,ibfd,i;
  for(i=0;i<nelem;i++)
    {
      iq=iptr[pc*i];
      ix=iptr[pc*i+1];
      ibv=iptr[pc*i+2];      
      ibvd=iptr[pc*i+3];
      ibf=iptr[pc*i+6];
      ibfd=iptr[pc*i+7];
      
      checkGradients(x+ix, q+iq,bv+ibv, bvd+ibvd,bf+ibf, bfd+ibfd,pde,d,e,p);
    }
}
     

void COMPUTE_RESIDUAL(double *R, double *mass, double *q, double *detJ, double *fflux,
		      double *bv, double *bvd,
		      double *bf, double *bfd,
		      int *iptr,int *iptrf,int *elem2face,
		      int pc, int pf,
		      int pde, int d, int e, int p, int nelem)
{
  int i,ix,idet,im,iR,iq,ibv,ibvd,ibf,ibfd;
  int nfp=facePerElem[e];
  for(i=0;i<nelem;i++)
    {
      //printf("elem %d:\n",i);
      ix=pc*i;
      iR=iq=iptr[ix];
      ibv=iptr[ix+2];
      ibvd=iptr[ix+3];
      idet=iptr[ix+5];
      ibf=iptr[ix+6];
      ibfd=iptr[ix+7];
      im=iptr[ix+10];
            
      volIntegral(R+iR,bv+ibv,bvd+ibvd,q+iq,detJ+idet, pde,d,e,p);
      faceIntegral(R+iR,fflux,bf+ibf,bfd+ibfd,elem2face+nfp*i,iptrf,q,pf,pde,d,e,p,i);
      invertMass(mass+im,R+iR,pde,d,e,p,i);
      //printf("%d %f\n",i,R[iR]);
    }
}


void COMPUTE_RHS(double *R,double *mass,double *bv, double *bvd, double *JinvV, double *detJ,
		 double *bf, double *bfd, double *JinvF,
		 double *faceWeight, double *fnorm, double *fflux,
		 double *x, double *q, int *elem2face, int *iptr, int *iptrf, int *faces,
		 int pc, int pf, int pde, int d , int e, int p, int nfaces, int nelem)
{

  FILL_FACES(fnorm, fflux, elem2face, iptr, iptrf, 
	     faceWeight, bf, bfd, q, 
	     pc, pf, pde, d, e, p, nelem,nfaces);

  FILL_BC(fnorm,fflux,faces,pde,d,e,p,nfaces);

  COMPUTE_FACE_FLUXES(fnorm,fflux,pde,d,e,p,nfaces,faces);

  //CHECK_GRADIENTS(x, q,bv, bvd, bf, bfd,iptr,pc,pde,d,e,p,nelem);

  COMPUTE_RESIDUAL(R,mass,q,detJ,fflux,
		   bv,bvd,bf,bfd,
		   iptr,iptrf,elem2face,
		   pc,pf,
		   pde,d, e, p, nelem);

}


void UPDATE_DOFS(double *qdest, double coef, double *qsrc, double *R, int ndof)
{
  int i;
  for(i=0;i<ndof;i++)
    qdest[i]=qsrc[i]+coef*R[i];
}
