double consIntegral(double *bv, double *q, double *detJ,
		    int pde, int d, int e, int p, int eid, int f1)
{
  int b,w,i,j,l,ld,m,f,n;
  int nbasis=order2basis[e][p];
  double wgt;
  double *bvv;
  int nfields=get_nfields[pde](d);
  double qv[nfields],qvv;
  int g=p2g[e][p];

  for(f=0;f<nfields;f++) qv[f]=0;
  
  for(w=0;w<ngElem[e][p];w++)
    {
      wgt=gauss[e][g][(d+1)*w+2]*detJ[w];
      bvv=bv+w*nbasis;
      for(f=0;f<nfields;f++)
	{
	  qvv=bvv[0]*q[f*nbasis];
	  for(b=1;b<nbasis;b++)     
	    qvv+=(bvv[b]*q[f*nbasis+b]);
	  qv[f]+=(wgt*qvv);
	}
    }
  return qv[f1];     
}


double consFaceIntegral( double *fflux,  int *elem2face,
			 int *iptrf,double *q, int pf, int pde, int d, int e, int p, int ielem,
			 int *faces, int eid, int f1)
{
  int b,w,i,j,l,ld,m,f,fid,fst,fsgn,neigh;
  int nbasis=order2basis[e][p];
  double wgt,v;
  double *bvv,*bvvd;
  int nfields=get_nfields[pde](d);
  double flux[nfields];
  int g=p2gf[e][p];
  int nfp=facePerElem[e];
  double resf[nfields];

  l=ld=m=0;
  for(f=0;f<nfields;f++)
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
	  wgt=gaussgl[e][g][(d)*w+1]*fsgn;
	  m=0;
	  if(neigh == -1){
	  for(f=0;f<nfields;f++)
	    {
	      // measure fluxes only at outer boundaries
	        //if (f==0) printf("flux=%f\n",flux[f]);
                flux[f]=fflux[fst+f]*wgt;
		resf[m]-=(flux[f]);
		m++;
	    }
	  } 
	  fst+=(3*fsgn*nfields);
	}
    }
  return resf[f1];
}


//XXX Debug here
double cutFaceCons(double *fflux, 
		   int pf, int pde, int d, int e, int p, int iorig,
		   int *cutoverset, int debug, int* cut2neigh, int* iblank, int f1)
{
  int b,w,i,j,m,f,floc,z;
  int nbasis=order2basis[e][p];
  double wgt,v;
  double *bvv,*bvvd;
  int nfields=get_nfields[pde](d);
  double flux;
  int g=p2gf[e][p];
  int nfp=facePerElem[e];
  int compute;
  double resf[nfields];
  double sgn; 

  for(f=0;f<nfields;f++)
    resf[f]=0;

  m=z=0;
  floc = 2*nfields; 
  for(i=0;i<nfp;i++)
    {
      //printf("\n"); 
      for(w=0;w<ngGL[e][p];w++)
	{
	  // subtract fluxes if it's the overset boundary
	  // add if it's a regular cut face
          if(cutoverset[z]>-1){
	    sgn = -1; 
	  }
	  else{
	    sgn = 1; 
	  }
	  wgt=gaussgl[e][g][(d)*w+1]*sgn; 
          // get the basis and basis derivative for this gauss point
	  m=0;
	  for(f=0;f<nfields;f++)
	    {
	      compute=(cut2neigh[i]==-1);
	      flux=fflux[floc+f]*wgt;
              //if (f==0) {
                //printf("cutface_flux : %d %d %e\n",i,cut2neigh[i],flux);
              //} 
	      //notice the sign change from the faceIntegral routine
	      if (compute) resf[m]+=(flux);
	      m++;
	    }
	  floc+=(3*nfields); 
	  z++; 
	}
    }
  //printf("iorig,resf[f1]=%d %.18e\n",iorig,resf[f1]);
  return resf[f1];
}

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
  int b,w,i,j,l,ld,m,f,fid,fst,fsgn,neigh,compute;
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
              compute=1;
              if (neigh >= 0) compute=(iblank[neigh]!=1);
	      if(compute){
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
	  printf("cut wgt = %f\n",wgt); 
	  printf("cut w = %i, qv = %f, qvd[%i] = %f %f\n",w,f,qv[f],qvd[f][0],qvd[f][1]);
	}
	for(f=0;f<nfields;f++){
	  printf("cut flux[%i] = %f %f\n",f,flux[0][f],flux[1][f]);
	}
	for(b=0;b<nbasis;b++){
	  printf("cut bvvd[%i] = %f %f\n",b,bvvd[0 + b*d], bvvd[1 + b*d]);
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
    for(f=0;f<nfields;f++)
      for(b=0;b<nbasis;b++)
	printf("R_volcut[%i] = %f\n",f*nbasis+b,keep[f*nbasis+b]);
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

//XXX Debug here
void cutFace(double *residual, double *fflux, double *bfL, double *bfR, 
	     double *q, int pf, int pde, int d, int e, int p, int iorig,
	     int *cutoverset, int debug, int* cut2neigh, int* iblank)
{
  int b,w,i,j,l,ld,m,f,floc,z;
  int nbasis=order2basis[e][p];
  double wgt,v;
  double *bvv,*bvvd;
  int nfields=get_nfields[pde](d);
  double flux;
  int g=p2gf[e][p];
  int nfp=facePerElem[e];
  int compute;
  double keep[nbasis*nfields],sgn; 

  for(f=0;f<nfields;f++)
    for(b=0;b<nbasis;b++)
      keep[f*nbasis+b]=0;

  l=ld=m=z=0;
  floc = 2*nfields; 
  for(i=0;i<nfp;i++)
    {
      //printf("\n"); 
      for(w=0;w<ngGL[e][p];w++)
	{
	  // subtract fluxes if it's the overset boundary
	  // add if it's a regular cut face
	  if(debug)	  printf("cutoverset = %i\n",cutoverset[z]);
          if(cutoverset[z]>-1){
	    sgn = -1; 
	  }
	  else{
	    sgn = 1; 
	  }
	  if(debug) printf("sgn = %f\n",sgn); 
	  wgt=gaussgl[e][g][(d)*w+1]*sgn; 
	  if(debug) printf("wgt = %f\n",wgt);
          // get the basis and basis derivative for this gauss point
          bvv=bfL+l;
          l+=nbasis;

	  m=0;
	  for(f=0;f<nfields;f++)
	    {
              compute=1;
              if (cut2neigh[i] >= 0) compute=(iblank[cut2neigh[i]]!=1);
              if(compute){
                flux=fflux[floc+f]*wgt;
		//if (cut2neigh[i]==-2 && f==0) printf("overset flux : %d %.18e\n",iorig,fflux[floc+f]);
  	        for(b=0;b<nbasis;b++)
	  	  {
 	            //notice the sign change from the faceIntegral routine
  		    residual[m]+=(flux*bvv[b]);
		    keep[m] += (flux*bvv[b]);
		    if(debug && f ==0)		  printf("side %i, w %i, field %i, basis %i,\n\tsgn = %f, flux = %f, bvv = %f, res inc = %f, curr res = %f\n",i,w,f,b,sgn,flux,bvv[b],flux*bvv[b],residual[m]);
		    m++;
		  }
	      }
	    }
	  floc+=(3*nfields); // increment through flux array, 3 is where the completed fluxes are computed
	  z++; 
	}
    }

  if(debug){
    printf("\n"); 
    m=0;
    for(f=0;f<nfields;f++){
      for(b=0;b<nbasis;b++){
	printf("R_facecut[%i] = %f\n",f*nbasis+b,keep[m]);
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
      else if (faces[6*i+4]==-2 ) {
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
			 int necut, int pccut, int *iptrc, double *fcflux, double *fwcut,int* cut2e,
			 int* cutoverset, int* cut2neigh,int* iblank)
{
  int nfields=get_nfields[pde](d);
  double normal[d],xnorm[d];
  int i,j,k,g,m,w,floc,wloc,ifl,ifr,iflux,ic2n,f,eid,ico,neigh;
  int nfp = facePerElem[e];

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
      */

      //if(i==1) printf("\tflx = %f %f %f %f\n\n",fflux[iflux+0],fflux[iflux+1],fflux[iflux+2],fflux[iflux+3]); 
    }
  //printf("\n======================\ngetting cut face fluxes\n============================\n");
  // cut face fluxes
  for(i=0;i<necut;i++){
    eid=cut2e[i];
    floc=iptrc[i*pccut+11];
    wloc=iptrc[i*pccut+9]; 
    ic2n=iptrc[i*pccut+12]; 
    ico=iptrc[i*pccut+13]; 
    int debug=0; 
    /*
    if(eid==916 ){
      debug = 1;
    }
    else{
      debug = 0;
    }
    */

    g = m = 0; 
    for(j=0;j<nfp;j++){
      for(w=0;w<ngGL[e][p];w++){

        for(k=0;k<d;k++){
	  if(cutoverset[ico+g]>-1){ // overset face needs sign flip
 	    xnorm[k]=-fwcut[wloc + m++];
          }
	  else{
 	    xnorm[k]=fwcut[wloc + m++];
          }

        }

        ifl=floc; 
        ifr=ifl+nfields;
        iflux=ifr+nfields;
	neigh = cut2neigh[ic2n+j];	
        debug=0;
	//if (neigh==-2) debug=1;
	if(debug){
	  printf("\nORIG %i, CUT i %i, j %i, w %i\n",eid,i,j,w);
	  printf("\tcutoverset = %i\n",cutoverset[ico+g]); 
	  printf("\tcut2neigh = %i\n",cut2neigh[ic2n+j]); 
	  printf("\tifl %i, ifr %i, iflux %i\n",ifl,ifr,iflux);
	  printf("\tLflx = %f %f %f %f\n",fcflux[ifl+0],fcflux[ifl+1],fcflux[ifl+2],fcflux[ifl+3]); 
	  printf("\tRflx = %f %f %f %f\n",fcflux[ifr+0],fcflux[ifr+1],fcflux[ifr+2],fcflux[ifr+3]); 
	  printf("\txNorm = %f %f \n",xnorm[0],xnorm[1]);
	}
	if(neigh!=eid && iblank[neigh]!=1) // ignore cut edges interior to orig elem
          gradient_indep_flux[pde](fcflux+ifl,fcflux+ifr,fcflux+iflux,xnorm,0.0);

	if(debug){
	  printf("\tflx%i = [%f %f %f %f]\n",(j*ngGL[e][p]+w),fcflux[iflux+0],fcflux[iflux+1],fcflux[iflux+2],fcflux[iflux+3]); 
	  printf("\tneigh = %i, eid = %i, blank = %i\n",neigh,eid,iblank[neigh]);
	}


	floc = floc + 3*nfields;
	g++; 
      }    // gauss pts
    } // tri faces
  } // cut elems

}

void invertMass(double *mass, double *R, int pde, int d , int e, int p,int iscut,int ireg, int debug)
{
  int i,j,f;
  int nbasis=order2basis[e][p];
  int iflag;
  int nfields=get_nfields[pde](d);

  if(iscut){
    if(ireg){
      if(debug) printf("\nEntering solvec_reg\n"); 

      solvec_copy_reshape_reg(mass,R,&iflag,nbasis,nfields,debug);  
      for(i=0;i<nbasis;i++)
	if(debug) printf("R_final[%i] = %f\n",i,R[i]);
    }
    else{
      if(debug) printf("\nEntering solvec\n"); 
      solvec_copy_reshape(mass,R,&iflag,nbasis,nfields);  
    }
  }
  else{
    if(debug) printf("\nEntering solvec\n"); 
    solvec_copy_reshape(mass,R,&iflag,nbasis,nfields);
  }
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
                      int *iptrc, int necut, int* cut2e, int* cut2neigh, int* iblank, int ireg, 
		      int *cutoverset,int imesh,int* faces)

{
  int i,j,k,ix,idet,im,iR,iq,ibv,ibvd,ibf,ibfd,eid,iflx,ic2n,ico;
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

      if(imesh==1 && (i==916) && 1==0){
	printf("DEBUG: Mesh %i, full cell %i:\n",imesh,i);
	debug = 1;
      }
      else if(imesh==0 && i==895 && 1==0){
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

      if(iblank[i]!=1){
        volIntegral(R+iR,bv+ibv,bvd+ibvd,q+iq,detJ+idet, pde,d,e,p,i);
	if(debug){
	  for(int f = 0; f<nfields; f++)
	    for(int j = 0; j<nbasis; j++)
	      printf("\tonly vol R(f = %i, b = %i) = %f\n",f,j,R[iR+f*nbasis+j]);
	}

        faceIntegral(R+iR,fflux,bf+ibf,bfd+ibfd,elem2face+nfp*i,iptrf,q,pf,pde,d,e,p,i,faces,iblank,i);
	if(debug){
	  for(int f = 0; f<nfields; f++)
	    for(int j = 0; j<nbasis; j++)
	      printf("\tfull R(f = %i, b = %i) = %f\n",f,j,R[iR+f*nbasis+j]);
	}

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
      ico=iptrc[ix+13];
      
      if(imesh==0 && eid==916 && 1==0){
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

      cutFace(R+iR,fcflux+iflx,bfcutL+ibf,bfcutR+ibf,q,pf,pde,d,e,p,eid,cutoverset+ico,debug,cut2neigh+ic2n,iblank);

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
      if(imesh==0 && (i==916||i==895) && (1==0)){
	debug = 1; 
	printf("MESH %i ELEM %i \n", imesh,i);
	for(int f = 0; f<nfields; f++)
	  for(int j = 0; j<nbasis; j++){
	    //R[iR+f*nbasis+j] = 0.001; 
	    printf("\tComplete R(f = %i, b = %i) = %f\n",f,j,R[iR+f*nbasis+j]);
	  }
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
	  break;    
        }
      invertMass(mass+im,R+iR,pde,d,e,p,iscut,ireg,debug);

      if(imesh==0 && (i==0||i==895) && (1==0)){
	printf("\n");
	for(int f = 0; f<nfields; f++)
	  for(int j = 0; j<nbasis; j++)
	    printf("\tUpdate(f = %i, b = %i) = %f\n",f,j,R[iR+f*nbasis+j]);
      }
    }
}

void COMPUTE_RHS(double *R,double *mass,double *bv, double *bvd, double *JinvV, double *detJ,
		 double *bf, double *bfd, double *JinvF,
		 double *faceWeight, double *fnorm, double *fflux,
		 double *x, double *q, int *elem2face, int *iptr, int *iptrf, int *faces,
		 int pc, int pf, int pccut, int pde, int d , int e, int p, int nfaces, int nelem,
                 double *bvcut, double *bvdcut,double *detJcut,
                 double *bfcutL, double *bfcutR,double *fwcut,
                 double *fcflux,int *iptrc,
                 int necut, int* cut2e, int *cut2face, int* cut2neigh, int* iblank, int ireg,
		 int* cutoverset, int imesh)
{

  FILL_FACES(x, fnorm, fflux, elem2face, iptr, iptrf, 
	     faceWeight, bf, bfd, q, 
	     pc, pf, pde, d, e, p, nelem,nfaces, 
	     fcflux, cut2face, cut2neigh, iptrc,  
             fwcut, bfcutL, bfcutR, cut2e, necut,pccut,faces,iblank);

  FILL_BC(fnorm,fflux,faces,pde,d,e,p,nfaces);

  COMPUTE_FACE_FLUXES(fnorm,fflux,pde,d,e,p,nfaces,faces,necut,pccut,iptrc,fcflux,fwcut,cut2e,cutoverset,cut2neigh,iblank);

  //CHECK_GRADIENTS(x, q,bv, bvd, bf, bfd,iptr,pc,pde,d,e,p,nelem);

  COMPUTE_RESIDUAL(R,mass,q,detJ,fflux,
		   bv,bvd,bf,bfd,
		   iptr,iptrf,elem2face,
		   pc,pf,pccut,
		   pde,d, e, p, nelem,
                   detJcut, fcflux, 
                   bvcut, bvdcut, bfcutL,  
		   bfcutR, 
                   iptrc, necut, cut2e, cut2neigh,iblank,ireg,cutoverset,imesh,faces); 


}


void UPDATE_DOFS(double *qdest, double coef, double *qsrc, double *R, int ndof)
{
  int i;
  for(i=0;i<ndof;i++)
    qdest[i]=qsrc[i]+coef*R[i];
}


double COMPUTE_CONSERVATION(double *q, double *detJ, double *bv, int *iptr, int *iblank,
			    int pc, int pde, int d, int e, int p, int fieldid,int nelem,
			    int *cut2e,int *iptrc, double *bvcut, double *detJcut, 
			    int pccut, int necut, double *fcflux, int *cutoverset, int *cut2neigh,
			    double *fflux, int *elem2face, int *faces, int *iptrf, int pf,double *faceFluxSum,double dt)
{
  int i,ix,iq,ibv,idet,iflx,ic2n,ico;
  int eid;
  int nfp=facePerElem[e];
  double cons=0;
  double fcons=0;
  double cons1;
  double ctot=0;
  for(i=0;i<nelem;i++)
    {
     if (iblank[i]!=1) {
      ix=pc*i;
      iq=iptr[ix];
      ibv=iptr[ix+2];
      idet=iptr[ix+5];
      cons+=consIntegral(bv+ibv,q+iq,detJ+idet,pde,d,e,p,i,fieldid);
      fcons+=(consFaceIntegral(fflux,elem2face+nfp*i,iptrf,q,pf,pde,d,e,p,i,faces,i,fieldid));
     }
    }
  //Modify removing contributions from cut cells
  double fcons1=0;
  for(i=0;i<necut;i++)
    {
      // get original element quantities
      eid = cut2e[i]; 
      //printf("eid=%d\n",eid);
      ix=pc*eid;
      iq=iptr[ix];	    
      //cut cell quantities
      ix=pccut*i; 
      ibv=iptrc[ix+2];
      idet=iptrc[ix+5];
      iflx=iptrc[ix+11];
      ic2n=iptrc[ix+12];
      ico=iptrc[ix+13];
      
      cons1=consIntegral(bvcut+ibv,q+iq,detJcut+idet,pde,d,e,p,eid,fieldid);      
      //printf("cons1=%.18e\n",cons1);
      fcons1+=cutFaceCons(fcflux+iflx,pf,pde,d,e,p,eid,cutoverset+ico,
      			 0,cut2neigh+ic2n,iblank,fieldid);
      cons-=cons1;
    }
  //printf("fcons1=%f %f %f\n",cons,fcons,fcons1);
  fcons+=fcons1;
  (*faceFluxSum)+=fcons;
  //printf("faceFluxSum=%f\n",*faceFluxSum);
  cons-=((*faceFluxSum)*dt);  
  //printf("cons=%f\n",cons);
  return cons;
}

