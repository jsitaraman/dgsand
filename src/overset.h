void setOversetFluxes(double *fcflux, double *bfcutR, double *qB, int* cut2neigh, int *iptrB, int d, int e, int p, int pc, int pde)
// computes R q values at cut cell interface
// using q from other mesh
{
  int nfp = facePerElem[e];
  int nbasis = order2basis[e][p];
  int ngauss = ngGL[e][p]; 
  int nfields=get_nfields[pde](d);
  int eid, qloc; 

  int floc = 0;
  int bloc = 0;  
  for(int j = 0; j<nfp; j++){
    // neg neigh id means it's on the other mesh
    eid = cut2neigh[j]; 
    for(int w=0; w<ngauss; w++){
      // notice the sign switch on the fluxes
      // compared to the other cut face fluxes. 
      // this flux is getting added to the system, 
      // not removed
      if(eid<0){
	qloc = iptrB[abs(eid)*pc]; 
	for(int f = 0; f<nfields; f++){
          fcflux[floc+f+nfields]=-bfcutR[bloc]*qB[qloc+f*nbasis];
          for(int b=1;b<nbasis;b++)
            fcflux[floc+f+nfields]-=bfcutR[bloc+b]*qB[qloc+f*nbasis+b];
        } // nfields
      } // eid
      bloc+=nbasis;
      floc+=3*nfields; // third set of values to be computed later, skip ahead to next quad pt
    } // ngauss
  } // nfp
}

void EXCHANGE_OVERSET(double* fcfluxA, double* bfcutRA, double* qB, int* iptrcA, int* iptrB, int* cut2neighA, int necutA, int pccut, int d, int e, int p, int pc, int pde)
{
  int ix, iq, ibf, ic2n, iflx;
  int nfp = facePerElem[e];
  int eid, flag; 

  // Loop over cut cells in mesh A
  for(int i = 0; i<necutA; i++){
    ix=pccut*i;
    iflx=iptrcA[ix+11]; 
    ic2n=iptrcA[ix+12];

    // does this cut cell have an overset boundary?
    flag = 0; 
    for(int j = 0; j<nfp; j++){
      // if neigh id<0, it's on mesh B
      eid = cut2neighA[ic2n+j];
      if(eid < 0) flag = 1; 
    } 
    if(flag){
      // cut cell quantities
      ibf=iptrcA[ix+6];

      // mesh B quantities
      iq=iptrB[pc*abs(eid)]; 

      // interpolate q fluxes from mesh B
      setOversetFluxes(fcfluxA+iflx, bfcutRA+ibf, qB+iq, cut2neighA+ic2n, iptrB, d, e, p, pc, pde);
    } // flag
  } // cut cells
}

void interpOversetCutNodes(double *xA, double *xB, int *iptrB, int pc,
			   int *cut2neighA, 
			   double *bfcutLA, double *bfcutRA, double *JinvB,
			   int d, int e, int p,int nelemB)


{
  int i,j,k,f,n,w,b,fid,eid,ixB,ind,ij;
  int nfp = facePerElem[e];
  int nbasis=order2basis[e][p];
  int ngauss = ngGL[e][p]; 
  double xloc,yloc,dx[2],u[2],rs[2], x0[2]; 

  int bloc = 0;
  for(j=0;j<nfp;j++){
    // find out if it's an overset boundary
    eid = cut2neighA[j];
    for(w=0;w<ngauss;w++){
      if(eid<0){

        // Find x coordinates of the desired face quadrature points
        // have bfcutR, use this to get xyz coords
        xloc = 0.0;
        yloc = 0.0;
        for(b=0;b<nbasis;b++){
          xloc = xloc + bfcutLA[b]*xA[b]; 
          yloc = yloc + bfcutLA[b]*xA[b+nbasis]; 
        }

	// loop over mesh B elements and find element that contains mesh A gauss pt
	int inside = 0;
	while(inside==0){
	  for(n=0; n<nelemB;n++){
	    ixB = iptrB[pc*n+1]; 

	    // get node 0 of mesh B element
	    u[0] = 0; 
	    u[1] = 0; 
	    for(int b=0; b<nbasis;k++){
	      x0[0]  += xB[b]       *basis[e][i](u);
	      x0[1]  += xB[b+nbasis]*basis[e][i](u);
	    } 
           
	    // use invJ of elem B to see if pt is inside element 
	    dx[0] = xloc-x0[0]; 
	    dx[1] = yloc-x0[1];
	    ij = iptrB[pc*n+4]; 
	    axb(JinvB+ij,dx,rs,2); 
           
	    // if inside, store mesh B elem number as negative
	    if(rs[0]>=0 && rs[0] <=1 && rs[1]>=0 && rs[1] <=1){
	      cut2neighA[j] = -n; 
	      inside = 1; 
	    }
	  } // mesh B elem

	  // if finished loop and still not found, something is wrong
          if((n+1)>nelemB){
            printf("ERROR! Can't find mesh B element for mesh A point (%f, %f)!\n",xloc,yloc); 
            exit(0); 
          }
        } // while not inside
      
	// get mesh B shape function values at quad pt and store in bfcutR
	for(int b=0; b<nbasis;k++)  bfcutRA[bloc+b] = basis[e][b](rs);

      } // if overset edge

      // move to next quad pt
      bloc+=nbasis;
    } // loop gauss
  } // loop edges
}


void SETUP_OVERSET(int* cut2neighA, int* iptrA, int* iptrB, int* iptrcA, int* iptrcB, double* xA, double* xB, double* xcutA, double* bfcutLA, double* bfcutRA, double* JinvB, int d, int e, int p, int pc, int pccut, int necutA, int nelemB)
// For overset boundaries in mesh A, find corresponding element and shape function values on mesh B
{
  int ix, iq, ibf, ibfd, ifw, ic2n,iflx, cibf,flag;
  int nfp = facePerElem[e];

  // Loop over cut cells in mesh A
  for(int i = 0; i<necutA; i++){
    flag = 0; 	   
    ix=pccut*i;
    ic2n=iptrcA[ix+12];

    // check to see if cut cell has overset boundary
    for(int j = 0; j<nfp; j++)  
      if(cut2neighA[ic2n+j] == -1) flag = 1; 

    if(flag){
      // cut cell quantities
      ibf=iptrcA[ix+6];
      ifw=iptrcA[ix+9];
      iflx=iptrcA[ix+11];
 
      // fill fcflux array on cut cell
      interpOversetCutNodes(xA+ix, xB, iptrB, pc, cut2neighA+ic2n, bfcutLA+ibf, bfcutRA+ibf, JinvB, d, e, p, nelemB);
    }
  }
}

