
void volIntegral(double *residual, double *bv, double *bvd, double *q, double *detJ,
		 int pde, int d, int e, int p, int eid)
{
  int b,w,i,j,l,ld,m,f,n;
  int nbasis=order2basis[e][p];
  double wgt;
  double *bvv,*bvvd;
  int nfields=get_nfields[pde](d);
  double flux[d][nfields];
  double qv[nfields],qvd[nfields][d],keep[nbasis*d*nfields];
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

/*if(eid==2){
printf("\n");
for(f=0;f<nfields;f++){
printf("wgt = %f\n",wgt); 
printf("w = %i, qv = %f, qvd[%i] = %f %f\n",w,f,qv[f],qvd[f][0],qvd[f][1]);
}
for(f=0;f<nfields;f++){
printf("flux[%i] = %f %f\n",f,flux[0][f],flux[1][f]);
}
for(b=0;b<nbasis;b++){
printf("bvvd[%i] = %f %f\n",b,bvvd[0 + b*d], bvvd[1 + b*d]);
}
}
*/

      m=0;
      n=0;
      for(f=0;f<nfields;f++)
        {
         for(j=0;j<d;j++) flux[j][f]*=wgt;
	 for(b=0;b<nbasis;b++)
	  {
	    if (w==0) residual[m]=0;
	    for(j=0;j<d;j++){
	      residual[m]+=(bvvd[b*d+j]*flux[j][f]); /* \grad b . F */
	      keep[n++]=(bvvd[b*d+j]*flux[j][f]);
	    }
	    m++;
	  }
        }

/*if(eid==2){
printf("\n");
m=0;
n=0; 
for(f=0;f<nfields;f++){
	 for(b=0;b<nbasis;b++){
	 for(j=0;j<d;j++){
	   printf("f %i, b %i, d %i, N %f, flux %f, vol res input = %f\n",f, b, j,bvvd[b*d+j], flux[j][f], keep[n++]);
}
}
printf("\n");
}
for(f=0;f<nfields;f++){
	 for(b=0;b<nbasis;b++)
		printf("w = %i, Rvol[%i] = %f\n",w,m,residual[m++]);
printf("\n");
}

}
*/
    }
}


void faceIntegral(double *residual, double *fflux, double *bf, double *bfd, int *elem2face,
		  int *iptrf,double *q, int pf, int pde, int d, int e, int p, int ielem,
		  int *faces, int* iblank, int eid)
{
  int b,w,i,j,l,ld,m,f,fid,fst,fsgn,neigh;
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
      neigh = faces[6*fid+2] == eid ? faces[6*fid+4] : faces[6*fid+2];
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
	      if(iblank[neigh]!=1){
                flux[f]=fflux[fst+f]*wgt;
  	        for(b=0;b<nbasis;b++)
	  	  {
                    resf[m]-=(flux[f]*bvv[b]);
  		    residual[m]-=(flux[f]*bvv[b]);
		    m++;
	 	  }
	        }
	    }

/*if(ielem==2){
printf("\n"); 
m = 0; 
  for(f=0;f<nfields;f++)
    for(b=0;b<nbasis;b++){
        printf("face %i, w %i, flux[%i] = %f, bvv[%i] = %f\n\tresf[%i]=%f\n",i, w,f,flux[f],b,bvv[b],m,flux[f]*bvv[b]); 
	m++;
     }
}
*/

	  fst+=(3*fsgn*nfields);
	}
    }
/*
if(ielem==2){
printf("\n"); 
m = 0; 
  for(f=0;f<nfields;f++)
    for(b=0;b<nbasis;b++){
        printf("f = %i, b = %i, face res =%f\n",f,b,resf[m]);
	m++;
     }
}
*/

}


void cutVol(double *residual, double *bv, double *bvd, double *q, double *detJ,
	    int pde, int d, int e, int p, int iorig, int debug)
{
  int b,w,i,j,l,ld,m,f;
  int nbasis=order2basis[e][p];
  double wgt;
  double *bvv,*bvvd;
  int nfields=get_nfields[pde](d);
  double flux[d][nfields];
  double qv[nfields],qvd[nfields][d];
  int g=p2g[e][p];
  double keep[nfields*nbasis];
  
  l=ld=0;

  for(f=0;f<nfields;f++)
  for(b=0;b<nbasis;b++)
  keep[f*nbasis+b]=0;

  // for all gauss-quadrature points 
  for(w=0;w<ngElem[e][p];w++)
    {
      wgt=gauss[e][g][(d+1)*w+2]*detJ[w];	  // is multiplying by detJcut correct?
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



 if(debug){
printf("\n");
for(f=0;f<nfields;f++){
printf("\tcut wgt = %f\n",wgt); 
printf("\tcut w = %i, qv = %f, qvd[%i] = %f %f\n",w,f,qv[f],qvd[f][0],qvd[f][1]);
}
for(f=0;f<nfields;f++){
printf("\tcut flux[%i] = %f %f\n",f,flux[0][f],flux[1][f]);
}
for(b=0;b<nbasis;b++){
printf("\tcut bvvd[%i] = %f %f\n",b,bvvd[0 + b*d], bvvd[1 + b*d]);
}
}

      for(f=0;f<nfields;f++)
        {
         for(j=0;j<d;j++) flux[j][f]*=wgt;
	 for(b=0;b<nbasis;b++)
	  {
	    for(j=0;j<d;j++){
 	      //notice the sign change from the volIntegral routine
	      residual[m]-=(bvvd[b*d+j]*flux[j][f]); /* \grad b . F */
              keep[m]+=(bvvd[b*d+j]*flux[j][f]);
	    }
	    m++;
	  }
        }
    }

if(debug){
printf("\n");
  for(f=0;f<nfields;f++){
  for(b=0;b<nbasis;b++)
  printf("\tR_volcut[%i] = %f\n",f*nbasis+b,keep[f*nbasis+b]);
  printf("\n");
  }
  }

}

  double signum(int a)
{
  if(a<0){
    return 1.0; 
  }
  else{
    return -1.0; 
  }
}

  void cutFace(double *residual, double *fflux, double *bfL, double *bfR, 
  	       double *q, int pf, int pde, int d, int e, int p, int iorig,
	       int *cutoverset, int debug, int* cut2neigh, int* cut2face, int* iblank,
	       double* OSFflux, double* OSFshpL, double* OSFxn, int OSFnseg)
{
  int b,w,i,j,l,ld,m,n,f,floc,z;
  int nbasis=order2basis[e][p];
  double wgt,v;
  double *bvv,*bvvd,*flux;
  int nfields=get_nfields[pde](d);
  int g=p2gf[e][p];
  int nfp=facePerElem[e];
  double keep[nbasis*nfields],sgn; 

  for(f=0;f<nfields;f++)
  for(b=0;b<nbasis;b++)
  keep[f*nbasis+b]=0;

  l=ld=m=z=0;
//  floc = 2*nfields; 
//
  for(i=0;i<nfp;i++)
    {

if(debug){  
m=0;
printf("\n");
            for(f=0;f<nfields;f++){
              for(b=0;b<nbasis;b++){
printf("\tstarting res[%i] = %f\n",m,residual[m]);
m++;
}
printf("\n");
}
}

if(debug)	  printf("\tcutoverset = %i, iblank neigh = %i, cut2face = %i\n",cutoverset[i],iblank[cut2neigh[i]],cut2face[i]);

      if(cutoverset[i]>-1){ // cut overset face
        for(n=0;n<OSFnseg;n++){
          for(w=0;w<ngGL[e][p];w++){
	    wgt=gaussgl[e][g][(d)*w+1]; 
            bvv = OSFshpL+n*ngGL[e][p]*nbasis+w*nbasis; 
            flux = OSFflux+n*ngGL[e][p]*3*nfields+w*3*nfields+2*nfields;

            m=0; 
            for(f=0;f<nfields;f++){
              for(b=0;b<nbasis;b++){
		keep[m]+=(wgt*flux[f]*bvv[b]);
                residual[m]-=(wgt*flux[f]*bvv[b]); // note minus sign b/c we're adding flux through overset face
if(f==0 && debug)		  printf("\tside %i, cut overset = 1\n\t seg %i, w = %f , field %i, basis %i, m = %i,\n\tOSFflux[%i] = %f, bvv[%i] = %f, res inc = %f, curr res = %f\n",i,n,wgt,f,b,m,n*ngGL[e][p]*3*nfields+w*3*nfields+2*nfields+f,flux[f],n*ngGL[e][p]*nbasis+w*nbasis+b,bvv[b],flux[f]*bvv[b]*wgt,residual[m]);
                m++;
	      }
            } // loop over fields

          } // loop over gauss
        } // loop over segs
        l += nbasis*ngGL[e][p]; // increment through basis aray
      }
      else if(cut2face[i]!=-1){ // skip if it's an internal cut face
      floc=i*ngGL[e][p]*3*nfields+2*nfields; // directly set floc
        for(w=0;w<ngGL[e][p];w++)
	  {
	    wgt=gaussgl[e][g][(d)*w+1]; 
            // get the basis and basis derivative for this gauss point
            bvv=bfL+l;
            l+=nbasis; 

            if(iblank[cut2neigh[i]]!=1){ // skip if neighbor is blanked
  	      m=0;
	      for(f=0;f<nfields;f++)
	        {
                  flux=fflux+floc;
  	          for(b=0;b<nbasis;b++)
	  	    {
		      keep[m]+=(wgt*flux[f]*bvv[b]);
  		      residual[m]+=(wgt*flux[f]*bvv[b]); // note plus sign b/c we're subtracting cut flux from entire edge flux
if(debug && f ==0)		  printf("\tside %i, cut overset = %i\n\tw %i = %f, field %i, basis %i, m %i,\n\tflux[%i] = %f, bvv = %f, res inc = %f, curr res = %f\n",i,cutoverset[i],w,wgt,f,b,m,floc+f,flux[f],bvv[b],flux[f]*bvv[b]*wgt,residual[m]);
		      m++;
		    }
	          }
	      }
	    floc+=(3*nfields); // increment through flux array, 3 is where the completed fluxes are computed
	    z++; 
	}
      } // if not cut face
    } // loop over faces


if(debug){
printf("\n"); 
m=0;
  for(f=0;f<nfields;f++){
  for(b=0;b<nbasis;b++){
  printf("\tR_facecut[%i] = %f\n",f*nbasis+b,keep[m]);
  m++ ;
  }
  printf("\n");
  }
}


}

void setFaceQuantities(double *fnorm,double *fflux,int *elem2face, int *iptrf,
		       double *faceWeight, double *bf,double *bfd, double *q, 
		       int nfields, int pf, int pde, int d , int e, int p, 
		       int* faces, int* iblank, int eid)
{
  int b,w,i,j,k,l,f,n,fid,floc,fst,fsgn,nst,neigh;
  int nbasis=order2basis[e][p];
  double qv[nfields];
  double *bvv;
  int nfp=facePerElem[e];
  l=k=0;

  for(i=0;i<nfp;i++)
    {
      fsgn=elem2face[i]/abs(elem2face[i]);
      fid=abs(elem2face[i])-1; 
      neigh = faces[6*fid+2] == eid ? faces[6*fid+4] : faces[6*fid+2];

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
	      if(iblank[neigh]!=1){
  	        fflux[floc]=bvv[0]*q[f*nbasis];
	        for(b=1;b<nbasis;b++)	    
	  	  fflux[floc]+=(bvv[b]*q[f*nbasis+b]); 
	      }
	      else{
  	        fflux[floc]=0.0;
	      }
	    }	  
	  for(j=0;j<d*fsgn;j++)
	    fnorm[nst+(n++)]=faceWeight[k+j];
          k+=d;
	  fst+=(3*fsgn*nfields); // incr forward if fsgn>0 and backward otherwise
	}
    }
}


void setCutFacesQuantities(double *x, double *q, int *iptr, int pc,
			   double *fcflux, 
			   int *cut2face, int *cut2neigh, double *fwcut, 
			   double *bfcutL, double *bfcutR,
			   int nfields, int pde, int d, int e, int p, 
			   int icut, int iorig, int* iblank)
// this routine accumulates L and R q values to subtract cut face fluxes 
{
  int i,j,k,f,w,b,fid,eid;
  int nfp = facePerElem[e];
  int nbasis=order2basis[e][p];
  int floc = 0;
  int bloc = 0; 

  for(j=0;j<nfp;j++){
    //need faceid of current uncut face
    fid=cut2face[j]; 
    for(w=0;w<ngGL[e][p];w++){
      if(iblank[cut2neigh[j]]!=1){
        // L side quantities
        eid = iorig;
        for(f=0;f<nfields;f++){ 
          fcflux[floc+f]=bfcutL[bloc]*q[iptr[eid*pc]+f*nbasis];
          for(b=1;b<nbasis;b++)
            fcflux[floc+f]+=bfcutL[bloc+b]*q[iptr[eid*pc]+f*nbasis+b];
        }// loop over nfields
    
        // R side quantities
        eid = cut2neigh[j];

        // only do if it's face cut by cut boundary interface
        // (overset boundaries already exchanged and internal cut faces cancel out)
        if(eid>-1 && fid!=-1){ 
          for(f=0;f<nfields;f++){ 
            fcflux[floc+f+nfields]=bfcutR[bloc]*q[iptr[eid*pc]+f*nbasis];
            for(b=1;b<nbasis;b++)
              fcflux[floc+f+nfields]+=bfcutR[bloc+b]*q[iptr[eid*pc]+f*nbasis+b];    
          }// loop over nfields
        } else if(eid==-1){ // inflow face
            far_field[pde](fcflux+floc+nfields);
        } // end of r side
      } 
      else{ // flux = 0 for faces with blanked neighbors
          // L side quantities
          for(f=0;f<nfields;f++){
            fcflux[floc+f]=0.0; // l side
  	    fcflux[floc+f+nfields]=0.0; // r side
	  }
      } // iblank
      bloc+=nbasis;
      floc+=3*nfields; // third set of values to be computed later, skip ahead to next quad pt

    } // loop over quad pts
  } // loop over faces
}

void FILL_FACES(double *x, double *fnorm, double *fflux, int *elem2face,int *iptr, int *iptrf,
		double *faceWeight, double *bf, double *bfd, double *q, 
		int pc, int pf, int pde, int d, int e, int p, int nelem, int nfaces,
		double *fcflux, int *cut2face, int *cut2neigh, int *iptrc,  
                double *fwcut, double *bfcutL, double *bfcutR, int* cut2e, int necut, int pccut,
		int* faces, int* iblank)
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
      setFaceQuantities(fnorm,fflux,elem2face+nfp*i, iptrf, 
			faceWeight+ifw, bf+ibf, bfd+ibfd, q+iq,
			nfields, pf, pde, d, e, p, faces, iblank,i);
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
      setCutFacesQuantities(x, q, iptr,pc, fcflux+iflx,
			    cut2face+ic2n,cut2neigh+ic2n,  
                            fwcut+ifw,bfcutL+ibf, bfcutR+ibf, 
			    nfields, pde, d, e, p, i, eid,iblank); 
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
			 int pde, int d, int e, int p, int nfaces, int *faces,
			 int necut, int pccut, int *iptrc, double *fcflux, double *fwcut,
                         int* OSFnseg, int* OSFeID, double* OSFflux, double* OSFxn, double* OSFshpL, double* OSFshpR, 
 			 int* cut2e, int* cutoverset, int* cut2neigh,int* iblank)
{
  int nfields=get_nfields[pde](d);
  double normal[d],xnorm[d];
  int i,j,k,g,m,n,w,floc,wloc,ifl,ifr,iflux,ic2n,f,eid,ico,neigh,ixn,iflx;
  int nfp = facePerElem[e];
  double *flux, *norm;
  double keep[3];

  // Full Faces
  for(i=0;i<nfaces*ngGL[e][p];i++)
    {
      f=i/ngGL[e][p];
      for(j=0;j<d;j++)
	    xnorm[j]=fnorm[d*i+j];

      ifl=i*3*nfields;
      ifr=ifl+nfields;
      iflux=ifr+nfields;
      gradient_indep_flux[pde](fflux+ifl,fflux+ifr,fflux+iflux,xnorm,0.0);

/*
if(i==1){
printf("FULL i %i, j %i, w %i\n",i,j,w);
printf("\tifl %i, ifr %i, iflux %i\n",ifl,ifr,iflux);
printf("\tLflx = %f %f %f %f\n",fflux[ifl+0],fflux[ifl+1],fflux[ifl+2],fflux[ifl+3]); 
printf("\tRflx = %f %f %f %f\n",fflux[ifr+0],fflux[ifr+1],fflux[ifr+2],fflux[ifr+3]); 
printf("\txNorm = %f %f \n",xnorm[0],xnorm[1]);
}


if(i==1) printf("\tflx = %f %f %f %f\n\n",fflux[iflux+0],fflux[iflux+1],fflux[iflux+2],fflux[iflux+3]); 
*/
    }
//printf("\n======================\ngetting cut face fluxes\n============================\n");

  // Cut Face fluxes
  int osfloc = 0;
  for(i=0;i<necut;i++){
    eid=cut2e[i];
    floc=iptrc[i*pccut+11];
    wloc=iptrc[i*pccut+9]; 
    ic2n=iptrc[i*pccut+12]; 
    ico=iptrc[i*pccut+12]; 
    ixn=iptrc[i*pccut+13]; 
    iflx=iptrc[i*pccut+16]; 

int debug; 
if(eid==1 && i==2){
debug = 1;
}
else{
debug = 0;
}
    g = m = 0; 

    for(j=0;j<nfp;j++){
      keep[0]=0.0; 
      keep[1]=0.0; 
      keep[2]=0.0; 
//if(debug) printf("mstart = %i\n",m);
      if(cutoverset[ico+j]>-1){ // cut overset face       

        for(n=0;n<OSFnseg[i];n++){
          for(w=0;w<ngGL[e][p];w++){
            for(k=0;k<d;k++) xnorm[k] = OSFxn[ixn + n*d*ngGL[e][p] + w*d + k]; // XXX double check index

	    flux=OSFflux + iflx + n*ngGL[e][p]*3*nfields + w*3*nfields;

            gradient_indep_flux[pde](flux,flux+nfields,flux+2*nfields,xnorm,0.0);
	    keep[0] = flux[0];
	    keep[1] = flux[0+nfields];
	    keep[2] = flux[0+2*nfields];
if(isnan(flux[2*nfields]) || isnan(flux[2*nfields+1]) ||isnan(flux[2*nfields+2]) ||isnan(flux[2*nfields+3])){
printf("\tflux index:\n\t\tiflux = %i\n\t\tseg ind = %i\n\t\tcur gauss = %i\n",iflx,n*ngGL[e][p]*3*nfields,w*3*nfields);
printf("\nORIG %i, CUT i %i, j %i, seg %i, w %i\n",eid,i,j,n,w);
printf("\tcutoverset = %i\n",cutoverset[ico+j]);
printf("\tcut2neigh = %i\n",cut2neigh[ic2n+j]);
printf("\tLflx = %f %f %f %f\n",flux[0],flux[1],flux[2],flux[3]); 
printf("\tRflx = %f %f %f %f\n",flux[4],flux[5],flux[6],flux[7]); 
//printf("\tLflx = %f \n",keep[0]);
//printf("\tRflx = %f \n",keep[1]);
printf("\txNorm[%i] = %f %f \n",ixn + n*d*ngGL[e][p] + w*d,xnorm[0],xnorm[1]);

printf("\tflx = %f \n",keep[2]);
exit(1);
}

if(debug){
printf("\tflux index:\n\t\tiflux = %i\n\t\tseg ind = %i\n\t\tcur gauss = %i\n",iflx,n*ngGL[e][p]*3*nfields,w*3*nfields);
printf("\nORIG %i, CUT i %i, j %i, seg %i, w %i\n",eid,i,j,n,w);
printf("\tcutoverset = %i\n",cutoverset[ico+j]); 
printf("\tcut2neigh = %i\n",cut2neigh[ic2n+j]); 
//printf("\tifl %i, ifr %i, iflux %i\n",ifl,ifr,iflux);
//printf("\tLflx = %f %f %f %f\n",fcflux[ifl+0],fcflux[ifl+1],fcflux[ifl+2],fcflux[ifl+3]); 
//printf("\tRflx = %f %f %f %f\n",fcflux[ifr+0],fcflux[ifr+1],fcflux[ifr+2],fcflux[ifr+3]); 
printf("\tLflx = %f \n",keep[0]);
printf("\tRflx = %f \n",keep[1]);
printf("\tflx = %f \n",keep[2]);
printf("\txNorm[%i] = %f %f \n",ixn + n*d*ngGL[e][p] + w*d,xnorm[0],xnorm[1]);
}



          }// loop over gauss pts
        }// loop over segments

	// increment counters for regular faces
	floc = floc + 3*nfields*ngGL[e][p]; 
        m = m+d*ngGL[e][p];
	g++;
      }
      else{ // regular cut face // XXX XNORM INDEX IS WRONG HERE
        for(k=0;k<d;k++) xnorm[k]=fwcut[wloc + m + k];
        //for(k=0;k<d;k++) xnorm[k]=fwcut[wloc + m++];
        for(w=0;w<ngGL[e][p];w++){
          ifl=floc; 
          ifr=ifl+nfields;
          iflux=ifr+nfields;

  	  neigh = cut2neigh[ic2n+j];
  	  if(neigh!=eid && iblank[neigh]!=1) // ignore cut edges interior to orig elem
            gradient_indep_flux[pde](fcflux+ifl,fcflux+ifr,fcflux+iflux,xnorm,0.0);
	    keep[0] = fcflux[ifl+2];
	    keep[1] = fcflux[ifr+2];
	    keep[2] = fcflux[iflux+2];

	  floc = floc + 3*nfields;
	  g++; 
/*
if(debug){
printf("\nORIG %i, CUT i %i, j %i, w %i\n",eid,i,j,w);
printf("\tcutoverset = %i\n",cutoverset[ico+j]); 
printf("\tcut2neigh = %i\n",cut2neigh[ic2n+j]); 
//printf("\tifl %i, ifr %i, iflux %i\n",ifl,ifr,iflux);
//printf("\tLflx = %f %f %f %f\n",fcflux[ifl+0],fcflux[ifl+1],fcflux[ifl+2],fcflux[ifl+3]); 
//printf("\tRflx = %f %f %f %f\n",fcflux[ifr+0],fcflux[ifr+1],fcflux[ifr+2],fcflux[ifr+3]); 
printf("\tLflx = %f \n",keep[0]);
printf("\tRflx = %f \n",keep[1]);
printf("\tflx = %f \n",keep[2]);
printf("\txNorm[%i] = %f %f \n",wloc,xnorm[0],xnorm[1]);
}
*/
        } // if gauss pt
	m = m+d*ngGL[e][p];
      } // if cutoverset


/*
if(debug){
printf("\tflx%i = [%f %f %f %f]\n",(j*ngGL[e][p]+w),fcflux[iflux+0],fcflux[iflux+1],fcflux[iflux+2],fcflux[iflux+3]); 
printf("\tneigh = %i, eid = %i, blank = %i\n",neigh,eid,iblank[neigh]);
}
*/
    } // tri faces
  } // cut elems

}


void invertMass(double *mass, double *R, int pde, int d , int e, int p,int iscut,int ireg, int debug, int ielem)
{
  int i,j,f;
  int nbasis=order2basis[e][p];
  int iflag;
  int nfields=get_nfields[pde](d);
  double b[nbasis];

  // store residual vector to check solution later
  for(int i=0;i<nbasis;i++) b[i] = R[i]; 

  if(iscut && ireg){

 printf("TEST2: R[ir] = %f\n",R[0]);

    solvec_copy_reshape_reg(mass,R,&iflag,nbasis,nfields,debug);  

if(debug){
 for(i=0;i<nbasis;i++)
 printf("update[%i] = %f\n",i,R[i]);
}
  }
  else{
    solvec_copy_reshape(mass,R,&iflag,nbasis,nfields);
  }

  // check accuracy of matrix solve
  checksol(mass,R,b,nbasis,ielem,debug);
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
                      int *iptrc, int necut, int* cut2e, int* cut2neigh, int* cut2face, int* iblank, int ireg, 
		      int *cutoverset,int imesh,int* faces,
                      double* OSFflux, double* OSFshpL, double* OSFxn, int* OSFnseg)
{
  int i,j,k,ix,idet,im,iR,iq,ibv,ibvd,ibf,ibfd,eid,iflx,ic2n,ixn,ishp,iflx2;
  int nfp=facePerElem[e];
  int f,w,ld,stop;
  double max;
  int iR2; 

int nfields=get_nfields[pde](d);
int nbasis=order2basis[e][p];
int debug;


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
 
iR2 = iptr[pc*(i+1)];
for(j=iR;j<iR2;j++)
if(isnan(R[j])){
printf("Pre-cut res for elem %i is NaN\n",i);
exit(1);
}
/*
if(imesh==1 && (i==916)){
printf("DEBUG: Mesh %i, full cell %i:\n",imesh,i);
debug = 1;
}
else if(imesh==0 && i==3){
printf("DEBUG: Mesh %i, full cell %i:\n",imesh,i);
debug = 1;
}
else{
debug = 0;
}
if(debug){ // print out node weights
for(int f = 0; f<nfields; f++)
for(int j = 0; j<nbasis; j++)
printf("\tq weights, q(f = %i, b = %i) = %f\n",f,j,q[iq+f*nbasis+j]);
}
*/

      if(iblank[i]!=1){
        volIntegral(R+iR,bv+ibv,bvd+ibvd,q+iq,detJ+idet, pde,d,e,p,i);
/*if(debug){
for(int f = 0; f<nfields; f++)
for(int j = 0; j<nbasis; j++)
printf("\tonly vol R(f = %i, b = %i) = %f\n",f,j,R[iR+f*nbasis+j]);
}
*/
        faceIntegral(R+iR,fflux,bf+ibf,bfd+ibfd,elem2face+nfp*i,iptrf,q,pf,pde,d,e,p,i,faces,iblank,i);
/*
if(debug){
for(int f = 0; f<nfields; f++)
for(int j = 0; j<nbasis; j++)
printf("\tfull R(f = %i, b = %i) = %f\n",f,j,R[iR+f*nbasis+j]);
}
*/

      }
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
      ic2n=iptrc[ix+12];
      ixn=iptrc[ix+13];
      ishp=iptrc[ix+15];
      iflx2=iptrc[ix+16];

if(imesh==1 && eid==1){
printf("DEBUG: Mesh %i, cut cell %i:\n",imesh,i);
printf("CUT VOL:\n"); 
debug = 1;
}
else{
debug = 0;
}
      cutVol(R+iR,bvcut+ibv,bvdcut+ibvd,q+iq,detJcut+idet,pde,d,e,p,eid,debug);

iR2 = iptr[pc*(eid+1)];
for(j=iR;j<iR2;j++)
if(isnan(R[j])){
printf("Post-vol-cut res for elem %i is NaN\n",eid);
exit(1);
}
if(debug) printf("\nCUT FACE: cut elem %i\n",i);

      cutFace(R+iR,fcflux+iflx,bfcutL+ibf,bfcutR+ibf,q,pf,pde,d,e,p,eid,cutoverset+ic2n,debug,cut2neigh+ic2n,cut2face+ic2n, iblank, OSFflux+iflx2, OSFshpL+ishp, OSFxn+ixn, OSFnseg[i]);

iR2 = iptr[pc*(eid+1)];
for(j=iR;j<iR2;j++)
if(isnan(R[j])){
printf("Post-face-cut res for elem %i is NaN\n",eid);
exit(1);
}


    }

  //Solve each element
max = 0; 
  for(i=0;i<nelem;i++)
    {
      ix=pc*i;
      iR=iq=iptr[ix];
      im=iptr[ix+10];

      stop = 0; 
      for(j=0;j<nbasis;j++)
	if(fabs(R[iR+j])>1e-10){
	  stop = 1;
	  if(fabs(R[iR+j]) > max) max = fabs(R[iR+j]);  
	}
#if 0      
      if(stop){
	printf("===============\n");
	printf("Debug elem %i\n",i);
	for(j=0;j<nbasis;j++) printf("R(%i) = %16.12e;\n",j+1,R[iR+j]); 
	printf("\n");
	ld = 0; 
	for(j=0;j<nbasis;j++) 
	  for(k=0;k<nbasis;k++){
	    printf("M(%i,%i) = %16.12e;\n",j+1,k+1, mass[im+ld]); 
	    ld++;
	  }
	
	printf("===============\n");
      }
#endif     

if(imesh==1 & (i==0)){
debug = 1; 
printf("MESH %i ELEM %i \n", imesh,i);
//for(int f = 0; f<nfields; f++)
f=0;
printf("TEST: R[ir] = %f\n",R[iR]);
for(int j = 0; j<nbasis; j++){
//R[iR+f*nbasis+j] = 0.001; 
printf("\tComplete R(f = %i, b = %i) = %f\n",f,j,R[iR+f*nbasis+j]);
}
printf("\n");
int m = 0; 
for(int k = 0; k<nbasis; k++)
for(int j = 0; j<nbasis; j++){
printf("mass(%i,%i) = %f\n",k+1,j+1,mass[im+m]);
m++;
}
}
else{
debug = 0; 
}

      int iscut = 0; 
      for(j=0;j<necut;j++)
        if(abs(i-cut2e[j])==0){
          iscut=1;
debug = 1; 
	  break;    
        }
//printf("Elem %i\n",i);
      invertMass(mass+im,R+iR,pde,d,e,p,iscut,ireg,debug,i);

//if(imesh==0 && (i==3||i==895)){
//printf("\n");
//for(int f = 0; f<nfields; f++)
//for(int j = 0; j<nbasis; j++)
//printf("\tUpdate(f = %i, b = %i) = %f\n",f,j,R[iR+f*nbasis+j]);
//}
//
    }
}

void COMPUTE_RHS(double *R,double *mass,double *bv, double *bvd, double *JinvV, double *detJ,
		 double *bf, double *bfd, double *JinvF,
		 double *faceWeight, double *fnorm, double *fflux,
		 double *x, double *q, int *elem2face, int *iptr, int *iptrf, int *faces,
		 int pc, int pf, int pccut, int pde, int d , int e, int p, int nfaces, int nelem,
                 double *bvcut, double *bvdcut,double *detJcut,
                 double *bfcutL, double *bfcutR,double *fwcut, double* fcflux,
                 int* OSFnseg, int* OSFeID, double* OSFxn, double* OSFshpL, double* OSFshpR, 
                 double* OSFflux, int *iptrc,
                 int necut, int* cut2e, int *cut2face, int* cut2neigh, int* iblank, int ireg,
		 int* cutoverset, int imesh)
{

  FILL_FACES(x, fnorm, fflux, elem2face, iptr, iptrf, 
	     faceWeight, bf, bfd, q, 
	     pc, pf, pde, d, e, p, nelem,nfaces, 
	     fcflux, cut2face, cut2neigh, iptrc,  
             fwcut, bfcutL, bfcutR, cut2e, necut,pccut,faces,iblank);

  FILL_BC(fnorm,fflux,faces,pde,d,e,p,nfaces);

  COMPUTE_FACE_FLUXES(fnorm,fflux,pde,d,e,p,nfaces,faces,necut,pccut,iptrc,fcflux,fwcut,
                      OSFnseg, OSFeID, OSFflux, OSFxn,OSFshpL, OSFshpR,
		      cut2e,cutoverset,cut2neigh,iblank);

  //CHECK_GRADIENTS(x, q,bv, bvd, bf, bfd,iptr,pc,pde,d,e,p,nelem);

  COMPUTE_RESIDUAL(R,mass,q,detJ,fflux,
		   bv,bvd,bf,bfd,
		   iptr,iptrf,elem2face,
		   pc,pf,pccut,
		   pde,d, e, p, nelem,
                   detJcut, fcflux, 
                   bvcut, bvdcut, bfcutL,  
		   bfcutR, 
                   iptrc, necut, cut2e, cut2neigh, cut2face,iblank,ireg,cutoverset,imesh,faces,
		   OSFflux, OSFshpL, OSFxn, OSFnseg); 
}

void UPDATE_DOFS(double *qdest, double coef, double *qsrc, double *R, int ndof)
{
  int i;
  for(i=0;i<ndof;i++)
    qdest[i]=qsrc[i]+coef*R[i];
}
