void setOversetFluxes(double *fcflux, double *bfcutR, double *qB, int* cutoverset, int *iptrB, int d, int e, int p, int pc, int pde, int debug)
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
if(debug) printf("In oversetfluxes\n");
  int m = 0; 
  for(int j = 0; j<nfp; j++){
    // neg neigh id means it's on the other mesh
    for(int w=0; w<ngauss; w++){
//printf("\t j = %i, w = %i\n",j,w);
      // notice the sign switch on the fluxes
      // compared to the other cut face fluxes. 
      // this flux is getting added to the system, 
      // not removed
      eid = cutoverset[m]; 
      if(eid>=0){


if(debug){
printf("\telemB = %i\n",eid); 
for(int f = 0; f<nfields; f++)
for(int k = 0; k<nbasis; k++)
printf("\tqB(f = %i, b = %i) = %f\n",f,j,qB[qloc+f*nbasis+k]);
}
	qloc = iptrB[eid*pc]; 
	for(int f = 0; f<nfields; f++){
          fcflux[floc+f+nfields]= bfcutR[bloc]*qB[qloc+f*nbasis];
          for(int b=1;b<nbasis;b++)
            fcflux[floc+f+nfields]+=bfcutR[bloc+b]*qB[qloc+f*nbasis+b];


/*if(isnan(bfcutR[bloc+b])){
  printf("\n\nbfcutR[%i] is NaN: %i %i %i %i \n\n",bloc+b,j,w,f,b);
exit(1); 
}
else if(isnan(qB[qloc+f*nbasis+b])){
  printf("\n\nqB[%i] is NaN: %i %i %i %i \n\n",qloc+f*nbasis+b,j,w,f,b);
exit(1); 
}
else{
printf("fcflux[%i] = %f\n",floc+f+nfields,fcflux[floc+f+nfields]);
}
}*/
//printf("\t\tfcflux[%i] = %f\n",floc+f+nfields,fcflux[floc+f+nfields]);

        } // nfields
      } // cut overset
      bloc+=nbasis;
      floc+=3*nfields; // third set of values to be computed later, skip ahead to next quad pt
      m++; 
    } // ngauss
  } // nfp
}

void EXCHANGE_OVERSET(double* fcfluxA, double* bfcutRA, double* qB, int* iptrcA, int* iptrB, int* cutoverset, int necutA, int pccut, int d, int e, int p, int pc, int pde, int imesh)
{
  int ip, ix, iq, ibf, ic2n, iflx, ico;
  int nfp = facePerElem[e];
  int eid, flag, debug, m;
//printf("In Exchange\n");
  // Loop over cut cells in mesh A
  for(int i = 0; i<necutA; i++){
    flag = 0; 
    ip=pccut*i;
    ix=iptrcA[ip+1];
    iflx=iptrcA[ip+11]; 
    ic2n=iptrcA[ip+12];
    ico = iptrcA[ip+13]; 

    // does this cut cell have an overset boundary?
    m = 0; 
    for(int j=0;j<nfp;j++)
      for(int w=0;w<ngGL[e][p];w++){
        if(cutoverset[ico+m]>=0){
  	eid = cutoverset[ico+m]; 
	  flag = 1;
        }
        m++;
      }

    if(flag ==1){
      // cut cell quantities
      ibf=iptrcA[ip+6];

      // mesh B quantities
      iq=iptrB[pc*eid]; 

if(imesh==1 && (i ==24 || i ==25) && (1==0)) {
debug=1;
}
else { debug = 0; }

      // interpolate q fluxes from mesh B
      setOversetFluxes(fcfluxA+iflx, bfcutRA+ibf, qB, cutoverset+ico, iptrB, d, e, p, pc, pde, debug);
    } // cut overset
  } // cut cells
}


double sign(double x1, double y1, double x2, double y2, double x3, double y3)
{
  return (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);
}
int pointInTri(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
// checks if point x0 y0 is in triangle defined by x and y 1-3
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
{
  double d1,d2,d3; 
  int neg,pos; 
  
  d1 = sign(x0,y0,x1,y1,x2,y2); 
  d2 = sign(x0,y0,x2,y2,x3,y3); 
  d3 = sign(x0,y0,x3,y3,x1,y1); 

  neg = (d1<0) || (d2<0) || (d3<0);
  pos = (d1>0) || (d2>0) || (d3>0);

  return !(neg && pos);
 
}
void interpOversetCutNodes(double *xA, double *xB, int *iptrB, int pc,
			   int *cutoversetA, 
			   double *bfcutLA, double *bfcutRA, double *JinvB,
			   int d, int e, int p,int nelemB,int debug)
{
  int i,j,k,f,n,w,b,fid,eid,ixB,ind,ij;
  int nfp = facePerElem[e];
  int nbasis=order2basis[e][p];
  int ngauss = ngGL[e][p]; 
  double xloc,yloc,dx[2],u[2],rs[2], x0[2]; 
  double xvert[6];

  int bloc = 0;
  int m; 
/*
if(debug){
// ================================================
printf("in interpOversetCutNodes\n"); 

//debug print vertices of cut cell
double xvert[6]; 
for(j=0;j<nfp;j++){
for(w=0;w<ngauss;w++){
xloc = 0; 
yloc = 0; 
	for(k=0;k<nbasis;k++){ // accumulate bases and vertex coords
          xloc   += bfcutLA[bloc+k]*xA[k];
          yloc   += bfcutLA[bloc+k]*xA[k+nbasis];
printf("\tbfcutL[%i] = %f \n",bloc+k,bfcutLA[bloc+k]);
        }
bloc = bloc+nbasis;
printf("face %i gauss %i= %f %f\n",j,w,xloc,yloc);
}
}

bloc = 0;

for(j=0;j<nfp;j++){
//print vertices
        u[0]=eloc[e][1][d*j];
        u[1]=eloc[e][1][d*j+1];
        xloc = 0.0;
        yloc = 0.0;
        for(k=0;k<nbasis;k++){ // accumulate bases and vertex coords
          xloc   += basis[e][k](u)*xA[k];
          yloc   += basis[e][k](u)*xA[k+nbasis];
printf("basis %i: %f\n",k,basis[e][k](u));
        }
printf("Vertex %i (%f,%f) = %f %f\n",j,u[0],u[1],xloc,yloc);
}
// ================================================
}

*/
  m = 0; 
  for(j=0;j<nfp;j++){
    for(w=0;w<ngauss;w++){
      if(cutoversetA[m]>=0){
        // Find x coordinates of the desired face quadrature points
        // have bfcutL, use this to get xyz coords
        xloc = 0.0;
        yloc = 0.0;
        for(b=0;b<nbasis;b++){
          xloc = xloc + bfcutLA[bloc+b]*xA[b]; 
          yloc = yloc + bfcutLA[bloc+b]*xA[b+nbasis]; 
        }

if(debug) printf("f %i, w %i, xloc,yloc = %f %f\n",j,w,xloc,yloc);

	// loop over mesh B elements and find element that contains mesh A gauss pt
	int inside = 0;
	  for(n=0; n<nelemB;n++){
	    ixB = iptrB[pc*n+1]; 

            // get 3 vertices of the current mesh B triangle
            for(int f=0;f<nfp;f++){
	      u[0]=eloc[e][1][d*f];
              u[1]=eloc[e][1][d*f+1];

	      xvert[2*f] = 0.0;
	      xvert[2*f+1] = 0.0;
              for(b=0;b<nbasis;b++){
		xvert[2*f]   +=  basis[e][b](u)*xB[ixB+b];
		xvert[2*f+1] +=  basis[e][b](u)*xB[ixB+b+nbasis];
              }
            }
            
	    inside = pointInTri(xloc,yloc,xvert[0],xvert[1],xvert[2],xvert[3],xvert[4],xvert[5]);	            
	    if(inside){
	      cutoversetA[m] = n; 
              break; 
	    }
	  } // mesh B elem

	  // compute the rst coordinate
	  dx[0] = xloc-xvert[0];
	  dx[1] = yloc-xvert[1];
          ij = iptrB[pc*cutoversetA[m]+4];
          axb(JinvB+ij,dx,rs,2);

	  // if finished loop and still not found, something is wrong
          if(inside==0 || rs[0]<0 || rs[0]>1 || rs[1]<0 || rs[1]>1){
            printf("ERROR! Can't find mesh B element for mesh A point (%f, %f)!\n",xloc,yloc); 
	    exit(0); 
          }
	  else{
if(debug)	    printf("\t found mesh A pt (%f, %f) in mesh B elem %i:  (%f, %f), (%f, %f), (%f, %f)\n",xloc,yloc,cutoversetA[m],xvert[0],xvert[1],xvert[2],xvert[3],xvert[4],xvert[5]);
if(debug)		printf("\t\t rst = %f %f\n",rs[0],rs[1]);
if(debug) 	    printf("\t cutoverset[%i] = %i\n",j,cutoversetA[m]);
          }
      
	// get mesh B shape function values at quad pt and store in bfcutR
	for(int b=0; b<nbasis;b++){
	  bfcutRA[bloc+b] = basis[e][b](rs);
	}

      } // if overset edge

      // move to next quad pt
      bloc+=nbasis;
      m++; 
    } // loop gauss
  } // loop edges
}


void SETUP_OVERSET(int* cut2e, int* cutoversetA, int* iptrA, int* iptrB, int* iptrcA, int* iptrcB, double* xA, double* xB, double* bfcutLA, double* bfcutRA, double* JinvB, int d, int e, int p, int pc, int pccut, int necutA, int nelemB)
// For overset boundaries in mesh A, find corresponding element and shape function values on mesh B
{
  int ip, ix, iq, ibf, ibfd, ifw, ic2n,iflx, cibf,flag, ico;
  int nfp = facePerElem[e];
  int eid, m; 

  // Loop over cut cells in mesh A
  for(int i = 0; i<necutA; i++){
    ip=pccut*i;
    eid = cut2e[i]; 
    ix  =iptrA[eid*pc+1];
    ibf =iptrcA[ip+6];
    ic2n=iptrcA[ip+12];
    ico=iptrcA[ip+13];
    flag = 0; 

    // check to see if cut cell has overset boundary
    m = 0; 
    for(int j = 0; j<nfp; j++)   
      for(int w = 0; w<ngGL[e][p]; w++){
        if(cutoversetA[ico+m]==1) flag = 1; 
        m++; 
      }

    if(flag==1){
      // cut cell quantities
 
      // fill fcflux array on cut cell
//      printf("ncut %i\n",i); 
int debug;       
if(eid==2 && (i==0 || i ==1) ){
printf("eid = %i, i = %i\n",eid,i);
debug = 1; 
} else{
debug = 0;}
      interpOversetCutNodes(xA+ix, xB, iptrB, pc, cutoversetA+ico, bfcutLA+ibf, bfcutRA+ibf, JinvB, d, e, p, nelemB,debug);

if(debug){
 printf("ico = %i\n", ico); 
 m = 0;
 for(int f=0;f<nfp;f++)
 for(int w=0;w<ngGL[e][p];w++){
   printf("\tdebug! f = %i, cutoverset = %i \n",f, cutoversetA[ico+m]);
   m++; 
 }
}



    }
  }
}

