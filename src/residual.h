
void volIntegral(double *residual, double *bv, double *bvd, double *q, double *detJ,
		 int pde, int d, int e, int p)
{
  int b,w,i,j,l,ld,m,f;
  int nbasis=order2basis[e][p];
  double wgt;
  double *bvv,*bvvd;
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
      bvv=bv+w*nbasis;
      bvvd=bvd+w*nbasis*d;
 
      /* project the q field to the gauss point */
      for(f=0;f<nfields;f++)
	{
	  qv[f]=bvv[0]*q[f*nbasis];
	  for(b=1;b<nbasis;b++)	    
	    qv[f]+=(bvv[b]*q[f*nbasis+b]);
	  for(j=0;j<d;j++)
	    {
	      qvd[f][j]=bvvd[j]*q[f*nbasis];
	      for(b=1;b<nbasis;b++)
		qvd[f][j]+=(bvvd[b*d+j]*q[f*nbasis+b]);
	    }
	}
      /* compute the flux function in each dimension */
      for(j=0;j<d;j++)
	flux_function[pde](flux[j],qv,qvd,j);
      m=0;
      for(f=0;f<nfields;f++)
        {
         for(j=0;j<d;j++) flux[j][f]*=wgt;
	 for(b=0;b<nbasis;b++)
	  {
	    if (w==0) residual[m]=0;
	    for(j=0;j<d;j++)
	      residual[m]+=(bvvd[b*d+j]*flux[j][f]); /* \grad b . F */
	    m++;
	  }
        }
    }
}


void faceIntegral(double *residual, double *fflux, double *bf, double *bfd, int *elem2face,
		  int *iptrf,double *q, int pf, int pde, int d, int e, int p, int ielem)
{
  int b,w,i,j,l,ld,m,f,fid,fst,fsgn;
  int nbasis=order2basis[e][p];
  double wgt,v;
  double *bvv,*bvvd;
  int nfields=get_nfields[pde](d);
  double flux[nfields];
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
      // make sure to get the right place to take the flux
      fst=iptrf[pf*(fid+(1-fsgn)/2)+1]-(1-fsgn)*nfields/2+(1+fsgn)*nfields;
      for(w=0;w<ngGL[e][p];w++)
	{
	  //v=gaussgl[e][g][(d)*w];
	  wgt=gaussgl[e][g][(d)*w+1]*fsgn;
          // get the basis and basis derivative for this gauss point
          bvv=bf+l;
          bvvd=bf+ld;
          l+=nbasis;
          ld+=(nbasis*d);

	  m=0;
	  for(f=0;f<nfields;f++)
	    {
              flux[f]=fflux[fst+f]*wgt;
	      for(b=0;b<nbasis;b++)
		{
                  resf[m]-=(flux[f]*bvv[b]);
		  residual[m]-=(flux[f]*bvv[b]);
		  m++;
		}
	    }
	  fst+=(3*fsgn*nfields);
	}
    }
}


void cutVol(double *residual, double *bv, double *bvd, double *q, double *detJ,
	    int pde, int d, int e, int p)
{
  int b,w,i,j,l,ld,m,f;
  int nbasis=order2basis[e][p];
  double wgt;
  double *bvv,*bvvd;
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
      bvv=bv+w*nbasis;
      bvvd=bvd+w*nbasis*d;
 
      /* project the q field to the gauss point */
      for(f=0;f<nfields;f++)
	{
	  qv[f]=bvv[0]*q[f*nbasis];
	  for(b=1;b<nbasis;b++)	    
	    qv[f]+=(bvv[b]*q[f*nbasis+b]);
	  for(j=0;j<d;j++)
	    {
	      qvd[f][j]=bvvd[j]*q[f*nbasis];
	      for(b=1;b<nbasis;b++)
		qvd[f][j]+=(bvvd[b*d+j]*q[f*nbasis+b]);
	    }
	}
      /* compute the flux function in each dimension */
      for(j=0;j<d;j++)
	flux_function[pde](flux[j],qv,qvd,j);
      m=0;
      for(f=0;f<nfields;f++)
        {
         for(j=0;j<d;j++) flux[j][f]*=wgt;
	 for(b=0;b<nbasis;b++)
	  {
	    for(j=0;j<d;j++)
 	      //notice the sign change from the volIntegral routine
	      residual[m]-=(bvvd[b*d+j]*flux[j][f]); /* \grad b . F */
	    m++;
	  }
        }
    }
}


//XXX Debug here
  void cutFace(double *residual, double *fflux, double *bfL, double *bfR, 
  	       double *q, int pf, int pde, int d, int e, int p, int ielem)
{
  int b,w,i,j,l,ld,m,f,floc;
  int nbasis=order2basis[e][p];
  double wgt,v;
  double *bvv,*bvvd;
  int nfields=get_nfields[pde](d);
  double flux;
  int g=p2gf[e][p];
  int nfp=facePerElem[e];

  l=ld=m=0;
  floc = 0; 
  for(i=0;i<nfp;i++)
    {
      for(w=0;w<ngGL[e][p];w++)
	{
	  wgt=gaussgl[e][g][(d)*w+1]; 
          // get the basis and basis derivative for this gauss point
          bvv=bfL+l;
          l+=nbasis;

	  m=0;
	  for(f=0;f<nfields;f++)
	    {
              flux=fflux[floc+f]*wgt;
	      for(b=0;b<nbasis;b++)
		{
 	          //notice the sign change from the faceIntegral routine
		  residual[m]+=(flux*bvv[b]);
		  m++;
		}
	    }
	  floc+=(3*nfields); // increment through flux array, 3 is where the completed fluxes are computed
	}
    }
}

void setFaceQuantities(double *fnorm,double *fflux,int *elem2face, int *iptrf,
		       double *faceWeight, double *bf,double *bfd, double *q, 
		       int nfields, int pf, int pde, int d , int e, int p)
{
  int b,w,i,j,k,l,f,n,fid,floc,fst,fsgn,nst;
  int nbasis=order2basis[e][p];
  double qv[nfields];
  double *bvv;
  int nfp=facePerElem[e];
  l=k=0;
for(i=0;i<nbasis;i++){
  printf("\tq = ");
  for(f=0;f<nfields;f++)
    printf("%f\t",q[f*nbasis+i]);
}
printf("\n");

  for(i=0;i<nfp;i++)
    {
      fsgn=elem2face[i]/abs(elem2face[i]);
      fid=abs(elem2face[i])-1; 
      // pick out the right location for inserting fields for this face
      // the cell with negative sign for the edge fills in backward order
      //
      // nst = fst if fsgn > 0 
      // fst = iptrf[pf*(fid+1)+1] - 2*nfields ? if fsgn < 0 
      nst=iptrf[pf*fid];
      fst=iptrf[pf*(fid+(1-fsgn)/2)+1]-(1-fsgn)*nfields;
      n=0;
      for(w=0;w<ngGL[e][p];w++)
	{
          bvv=bf+l;
          l+=nbasis;
	  for(f=0;f<nfields;f++)
	    {
	      floc=fst+f; // increments through fields and gauss pts
	      fflux[floc]=bvv[0]*q[f*nbasis];
	      for(b=1;b<nbasis;b++)	    
		fflux[floc]+=(bvv[b]*q[f*nbasis+b]); 
	    }	  
	  for(j=0;j<d*fsgn;j++)
	    fnorm[nst+(n++)]=faceWeight[k+j];
          k+=d;
	  fst+=(3*fsgn*nfields); // incr forward if fsgn>0 and backward otherwise
	}
    }
}


void setCutFacesQuantities(double *x, double *q, int *iptr, int pc,
			   double *fcnorm, double *fcflux, 
			   int *cut2face, int *cut2neigh, double *fwcut, 
			   double *bfcutL, double *bfcutR,
			   int nfields, int pde, int d, int e, int p, int icut, int iorig)
// this routine accumulates L and R q values to subtract cut face fluxes 
{
  int i,j,f,w,b,fid,eid;
  int nfp = facePerElem[e];
  int nbasis=order2basis[e][p];
  int floc = 0; 
  for(j=0;j<nfp;j++){
    //need faceid of current uncut face
    fid=cut2face[j]; 
    for(w=0;w<ngGL[e][p];w++){
      printf("\n\tface %i, gauss %i\n",j,w);
      printf("\t\tqL[%i-%i] = ",floc,floc+nfields-1);

      // L side quantities
      eid = iorig;
      for(f=0;f<nfields;f++){ 
        fcflux[floc+f]=bfcutL[0]*q[iptr[eid*pc]+f*nbasis];
        for(b=1;b<nbasis;b++)
          fcflux[floc+f]+=bfcutL[b]*q[iptr[eid*pc]+f*nbasis+b];
printf("%f\t",fcflux[floc+f]);
      }// loop over nfields
printf("\n");
 
    
      // R side quantities
      eid = cut2neigh[j];
printf("R neigh = orig elem %i\n",eid);
for(b=0;b<nbasis;b++)
printf("bR[%i] = %f\n",b,bfcutR[b]);

      printf("\t\tqR[%i-%i] = ",floc+nfields,floc+2*nfields);
      if(eid!=-1 && fid!=-1){ // if it's an internal cut face
        for(f=0;f<nfields;f++){ 
          fcflux[floc+f+nfields]=bfcutR[0]*q[iptr[eid*pc]+f*nbasis];
          for(b=1;b<nbasis;b++)
            fcflux[floc+f+nfields]+=bfcutR[b]*q[iptr[eid*pc]+f*nbasis+b];    
printf("%f\t",fcflux[floc+f+nfields]);
        }// loop over nfields
      }
      else{ // force overset fluxes to be inflow (will replace later)
	  far_field[pde](fcflux+floc+nfields);    
for(f=0;f<nfields;f++)
  printf("%f\t",fcflux[floc+f+nfields]);
      } // end of r side
printf("\n");

      floc+=3*nfields; // third set of values to be computed later, skip ahead to next quad pt

    } // loop over quad pts
  } // loop over faces



}

void FILL_FACES(double *x, double *fnorm, double *fflux, int *elem2face,int *iptr, int *iptrf,
		double *faceWeight, double *bf, double *bfd, double *q, 
		int pc, int pf, int pde, int d, int e, int p, int nelem, int nfaces,
		double* fcnorm, double *fcflux, int *cut2face, int *cut2neigh, int *iptrc,  
                double *fwcut, double *bfcutL, double *bfcutR, int* cut2e, int necut, int pccut)
{
  int i,ix,eid,ic2n;
  int ifw,ibf,ibfd,iq,iflx;
  int nfields=get_nfields[pde](d);
  int nfp=facePerElem[e];
  for(i=0;i<nelem;i++)
    {
      ix=pc*i;
      iq=iptr[ix]; 
      ibf=iptr[ix+6];
      ibfd=iptr[ix+7];
      ifw=iptr[ix+9];

      // have every elem fill in it's L state
      // do q = sum(N*q)      
      printf("\tOrig elem %i\n",i);
      setFaceQuantities(fnorm,fflux,elem2face+nfp*i, iptrf, 
			faceWeight+ifw, bf+ibf, bfd+ibfd, q+iq,
			nfields, pf, pde, d, e, p);
    }
  if(necut>0){
    for(i=0;i<necut;i++){
      eid=cut2e[i]; 

      // cut cell quantities
      ix=pccut*i;
      ibf=iptrc[ix+6]; 
      ifw=iptrc[ix+9]; 
      iflx=iptrc[ix+11];
      ic2n=iptrc[ix+12]; 

      // grab both L and R fluxes for cut cells
      printf("\tCut elem %i\n",i);
      setCutFacesQuantities(x, q, iptr,pc,fcnorm, fcflux+iflx,
			    cut2face+ic2n,cut2neigh+ic2n,  
                            fwcut+ifw,bfcutL+ibf, bfcutR+ibf, 
			    nfields, pde, d, e, p, i, eid); 
    }
  }
}

void FILL_BC(double *fnorm,double *fflux, int *faces,
	     int pde,int d, int e, int p, int nfaces)
{
  int i,j,ifl,ifp,inr;
  int nfields=get_nfields[pde](d);
  for(i=0;i<nfaces;i++)
    {
      /* TODO                                        */
      /* there should be a way to avoid this if loop */
      /* collect boundary faces in to a separate list in preproc */ 
      if (faces[6*i+4] == -1) {
	    for(j=0;j<ngGL[e][p];j++)
	      {
	      ifp=(i*ngGL[e][p]+j)*3*nfields+nfields;
	      far_field[pde](fflux+ifp);
	      }
      }
      else if (faces[6*i+4]==-2) {
          for(j=0;j<ngGL[e][p];j++)
            {
              inr=((i*ngGL[e][p]+j)*d);
              ifl=(i*ngGL[e][p]+j)*3*nfields;
              ifp=ifl+nfields;
              wall_bc[pde](fflux+ifp,fflux+ifl,fnorm+inr,d);
              //far_field[pde](fflux+ifp);
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

  // XXX need to add cut cells

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

  solvec_copy_reshape(mass,R,&iflag,nbasis,nfields);
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
	}
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
		      int pc, int pf, int pccut,
		      int pde, int d, int e, int p, int nelem,
                      double *detJcut, double *fcflux,
                      double *bvcut, double *bvdcut, 
		      double *bfcutL, double *bfcutR,
                      int *iptrc, int necut, int* cut2e, int* cut2neigh)

{
  int i,ix,idet,im,iR,iq,ibv,ibvd,ibf,ibfd,eid,iflx,ic2n;
  int nfp=facePerElem[e];

  //Accumulate residual on all elements
  for(i=0;i<nelem;i++)
    {
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
     // invertMass(mass+im,R+iR,pde,d,e,p,i);
    }

  //Modify residual for cut cells
  for(i=0;i<necut;i++)
    {
      // get original element quantities
      eid = cut2e[i]; 
      ix=pc*eid;
      iR=iq=iptr[ix];

      //cut cell quantities
      ix=pccut*i; 
      ibv=iptrc[ix+2];
      ibvd=iptrc[ix+3];
      idet=iptrc[ix+5];
      ibf=iptrc[ix+6];
      ibfd=iptrc[ix+7];
      iflx=iptrc[ix+11];
      
      cutVol(R+iR,bvcut+ibv,bvdcut+ibvd,q+iq,detJcut+idet,pde,d,e,p);
      cutFace(R+iR,fcflux+iflx,bfcutL+ibf,bfcutR+ibf,q,pf,pde,d,e,p,i);
    }

  //Solve each element
  for(i=0;i<nelem;i++)
    {
      ix=pc*i;
      iR=iq=iptr[ix];
      im=iptr[ix+10];
      invertMass(mass+im,R+iR,pde,d,e,p,i);
    }
}

void COMPUTE_RHS(double *R,double *mass,double *bv, double *bvd, double *JinvV, double *detJ,
		 double *bf, double *bfd, double *JinvF,
		 double *faceWeight, double *fnorm, double *fflux,
		 double *x, double *q, int *elem2face, int *iptr, int *iptrf, int *faces,
		 int pc, int pf, int pccut, int pde, int d , int e, int p, int nfaces, int nelem,
                 double *bvcut, double *bvdcut,double *detJcut,
                 double *bfcutL, double *bfcutR,double *fwcut,
                 double *fcnorm,double *fcflux,double *xcut, int *iptrc,
                 int necut, int* cut2e, int *cut2face, int* cut2neigh)


{

  FILL_FACES(x, fnorm, fflux, elem2face, iptr, iptrf, 
	     faceWeight, bf, bfd, q, 
	     pc, pf, pde, d, e, p, nelem,nfaces, 
	     fcnorm, fcflux, cut2face, cut2neigh, iptrc,  
             fwcut, bfcutL, bfcutR, cut2e, necut,pccut);

  FILL_BC(fnorm,fflux,faces,pde,d,e,p,nfaces);

  COMPUTE_FACE_FLUXES(fnorm,fflux,pde,d,e,p,nfaces,faces);

  //CHECK_GRADIENTS(x, q,bv, bvd, bf, bfd,iptr,pc,pde,d,e,p,nelem);

  COMPUTE_RESIDUAL(R,mass,q,detJ,fflux,
		   bv,bvd,bf,bfd,
		   iptr,iptrf,elem2face,
		   pc,pf,pccut,
		   pde,d, e, p, nelem,
                   detJcut, fcflux, 
                   bvcut, bvdcut, bfcutL,  
		   bfcutR, 
                   iptrc, necut, cut2e, cut2neigh); 


}


void UPDATE_DOFS(double *qdest, double coef, double *qsrc, double *R, int ndof)
{
  int i;
  for(i=0;i<ndof;i++)
    qdest[i]=qsrc[i]+coef*R[i];
}
