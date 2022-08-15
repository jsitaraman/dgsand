#include<stdio.h>
#include<stdlib.h>
#include<math.h>
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

if(imesh==1 && (i ==24 || i ==25)) {
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

/*
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
*/
int isBetween(double xA, double yA, double xB, double yB, double xC, double yC)
// checks to see if (xC,yC) is on the straight line between (xA,yA) and (xB,yB)
{
  double d1 = sqrt((xA-xC)*(xA-xC) + (yA-yC)*(yA-yC));
  double d2 = sqrt((xB-xC)*(xB-xC) + (yB-yC)*(yB-yC));
  double d3 = sqrt((xA-xB)*(xA-xB) + (yA-yB)*(yA-yB));

  double tol = 1e-12; 
  if(fabs(d3-d1-d2)<=tol){
    return 1;
  }
  else{
    return 0; 
  }

}

int isSame(double A, double B){
  double tol = 1e-14; 
  if(fabs(A-B)<tol){
    return 1;
  } else{
    return 0; 
  }
}

void orderCoords(double* in, double* out, int n)
// brute force way of reordering things
// Life would be easier if I could just use vectors in here
{
// in = initial array (2*n x 1)
// out = outal output array 
// n = number of pts


//debug
//printf("==========\nInside orderCords\n\n");
//printf("n = %i,\nIn coords:\n",n);
//
//This line is causing/stopping segfault somehow
//for(int i=0;i<n;i++) printf("\t%f %f\n",in[2*i],in[2*i+1]);

  // Assuming first 2 entries into in array are x and y of one of the end pts
  // Get parametric definition, [x y] = [x0 y0] + [nx ny]*s, s = 0-1
  double nx = in[2]-in[0];
  double ny = in[3]-in[1];
  double L = sqrt(nx*nx+ny*ny);

  out[0] = in[0];
  out[1] = in[1];
  out[2*n] = in[2];
  out[2*n+1] = in[3];
  double s1,s2,ds,mindiff,x,y;
  int ind;
 
  //loop over out
  for(int i=0;i<n-1;i++){
     y = out[2*i+1];
     s1 = (y-out[1])/ny;
     
     // loop over in
     mindiff = 1.0; 
     for(int j=0;j<n;j++){
       y = in[2*j+1];
       s2 = (y-out[1])/ny;
       ds = s2-s1;
//debug
//printf("i = %i, j = %i: s1 = %f, s2 = %f\n",i,j,s1,s2); 

 //      printf("fabs(ds) = %f, ds = %f\n",fabs(ds),ds);
       if(fabs(ds)<mindiff && ds > 0){ 
         ind = j;
       }
   //    printf("\tind = %i\n",ind); 
     } // loop over in

     //store entry in out
     out[2*(i+1)]   = in[2*ind];
     out[2*(i+1)+1] = in[2*ind+1];

  } // loop over out    
//debug
//printf("n = %i,\nOrdered coords:\n",n);
//for(int i=0;i<n;i++) printf("\t%f %f\n",out[2*i],out[2*i+1]);
//printf("\nExiting orderCords\n============\n");
}

void createOversetSegs(double* xseg, double* JinvA, int necutB, double* xcutB, int* cut2eB, int* iptrcB, int* OSFnseg, double* OSFwgt, double* OSFshpL, int pccut,int e, int p, int d, int debug)
// Break up overset face into segments, one for each overlapping element, and distribute quad pts
{
  int ip, eid, ixc; 
  int nfp = facePerElem[e];
  int i,j,k,m,uniq;
  int nbasis=order_to_basis(e,p);        // basis for solution
  int maxseg = 20; 
  double u[2], bv[nbasis],pt[2],pt2[2];
  int ngGL=get_ngGL(e,p);
  double xcpy[maxseg*ngGL*d]; // tmp copy of xseg

  // copy overset end pt coords, rest should be blank anyways
  for(i=0;i<2*d;i++) xcpy[i]=xseg[i]; 

  // loop over cut cells on mesh B
  for(i=0;i<necutB;i++){        
    ip=pccut*i;
    eid = cut2eB[i];
    ixc = iptrcB[ip+1];
    
    // Get vertices of cut triangle B
    for(j=0;j<nfp;j++){ // loop over vertices
      u[0]=eloc[e][1][d*j]; 
      u[1]=eloc[e][1][d*j+1];
 
      pt[0] = 0;
      pt[1] = 0;
      for(k=0;k<nbasis;k++){

        bv[k] = basis[e][k](u);  
        pt[0] += bv[k]*xcutB[ixc+k];
        pt[1] += bv[k]*xcutB[ixc+k+3];
      }
      // if vert is on cut overset, then store mesh B original cell ID and vertex coordinate
      if(isBetween(xcpy[0],xcpy[1],xcpy[2],xcpy[3],pt[0],pt[1])){

        // store vertex coordinate if unique
        uniq = 1;
        for(k=0;k<(OSFnseg[0]+1);k++){          
          if(isSame(xcpy[2*k],pt[0]) && isSame(xcpy[2*k+1],pt[1])) uniq = 0;
        }
	
        if(uniq){
          OSFnseg[0]++;                   	  	  
          xcpy[2*OSFnseg[0]]   = pt[0];
          xcpy[2*OSFnseg[0]+1] = pt[1];
        }
      } // found vertex in between
    } // loop over edges B
  } // loop over cut cell B

  // Reorder list of segment end pts if necessary
  if(OSFnseg[0]>1){
    orderCoords(xcpy,xseg,OSFnseg[0]+1);
  } else{
    xseg=xcpy; 
  }

// debug
if(debug){
printf("\n\t%i segments:\n",OSFnseg[0]);
for(i=0;i<OSFnseg[0]+1;i++)
printf("\t\t%f %f\n",xseg[2*i],xseg[2*i+1]);
}


  // Distribute gauss pts between each of the segments
  double nx, ny, s, L;
  m = 0; 
  for(i=0;i<OSFnseg[0];i++){
    // define segment end pts
    pt[0] = xseg[2*i];   
    pt[1] = xseg[2*i+1];   
    pt2[0] = xseg[2*(i+1)];   
    pt2[1] = xseg[2*(i+1)+1];   
    nx = pt2[0]-pt[0];
    ny = pt2[1]-pt[1];
    L = sqrt(nx*nx+ny*ny); 

    // distribute gauss pts
    int g=p2gf[e][p];
    for(j=0;j<ngGL;j++){
      s = gaussgl[e][g][2*j]; 

// don't need xgseg if i have OSFshpL. Can always recompute them
//      xgseg[2*(i*ngGL + j)]   = pt[0] + nx*s;		// gauss pt xy coord
//      xgseg[2*(i*ngGL + j)+1] = pt[1] + ny*s;

      // reusing pt2 for each of the gauss pts
      pt2[0] = pt[0] + nx*s;		// gauss pt xy coord
      pt2[1] = pt[1] + ny*s;
      OSFwgt[i*ngGL + j] = L*gaussgl[e][g][j*2+1]; 	// gauss pt weight

      // convert physical coord into full element rst coord
      // rst = [JinvA][xy]
      axb(JinvA,pt2,u,d);
      for(k=0;k<nbasis;k++){
	OSFshpL[m] = basis[e][k](u);
        m++;
      }

      

//printf("\t\tgauss pt %i at [%f %f] weight %f\n",i*ngGL+j,xgseg[2*(i*ngGL + j)],xgseg[2*(i*ngGL + j)+1],OSFwgt[i*ngGL + j]);
printf("\t\tgauss pt %i at [%f %f] weight %f\n",i*ngGL+j,pt2[0],pt2[1],OSFwgt[i*ngGL + j]);
    }
  }  // loop over overset segments
}

/*
void interpOverSegs()
// Find node number and shape function value across all the shape functions on overset segment
{

}
*/
// NEW ONE 
void SETUP_OVERSET(int* cut2e, int* cut2eB, int* cutoversetA, int* iptrA, int* iptrB, int* iptrcA, int* iptrcB, double* xA, double* xB, double* xcutA, double* xcutB, double* bfcutLA, double* bfcutRA, double* JinvA, double* JinvB, int* OSFnseg, int* OSFeID, double* OSFwgt, double* OSFshpL, double* OSFshpR, int d, int e, int p, int pc, int pccut, int necutA, int necutB)
// For overset boundaries in mesh A, find corresponding elements and shape function values on mesh B. Build multiple segments along overset boundary to handle discontinuous overset flux across multiple elements
{
  int ip, ix, ixc, iq, ibf, ibfd, ifw, ic2n,iflx, cibf,flag, ico, cid, iosf, ishp;
  int nfp = facePerElem[e];
  int eid, i,j,k,m, count;
  int nbasis=order_to_basis(e,p);        // basis for solution
  int ngGL=get_ngGL(e,p);
  double pt[2],xvertA[6],u[2],bv[nbasis];
  int maxseg=20; 
  double *xseg;
//  double *xgseg;
  xseg=dgsand_alloc(double,((maxseg+1)*d)); // physical coordinates of segment end pts
//  xgseg=dgsand_alloc(double,(maxseg*ngGL*d)); // physical coordinates of segment quad pts

  // Loop over cut cells in mesh A
  for(int i = 0; i<necutA; i++){
    ip=pccut*i;
    eid = cut2e[i]; 
    ix  =iptrA[eid*pc+1];
    ixc = iptrcA[ip+1]; 
    ibf =iptrcA[ip+6];
    ic2n=iptrcA[ip+12];
    ico=iptrcA[ip+13];
    iosf=iptrcA[ip+14];
    ishp=iptrcA[ip+15];
    flag = 0; 

    // only handle elements with overset boundaries
    m = 0;     
    flag = 0; 
    for(int j = 0; j<nfp; j++)   
      for(int w = 0; w<ngGL; w++){
        if(cutoversetA[ico+m]==1){
          flag = 1; 
          cid = m; 
	}
        m++; 
      }
    if(flag){

int debug;       
//if(eid==2 && (i==0 || i ==1) ){
printf("\neid = %i, cut id = %i, overset gauss= %i\n",eid,i,cid); 
debug = 1; 
//} else{
//debug = 0;}

      // reset counters and arrays for new segment
      for(k=0;k<(maxseg+1)*d;k++) xseg[k] = 0.0;
      OSFnseg[i] = 1; 

      // Get coord of gauss pt on overset boundary
      //
      pt[0] = 0.0; 
      pt[1] = 0.0;       
      for(int k=0;k<nbasis;k++){
        pt[0] += bfcutLA[ibf + cid*nbasis + k]*xA[ix+k];
        pt[1] += bfcutLA[ibf + cid*nbasis + k]*xA[ix+k+3];
      }

      // Get vertices of cut triangle A
      for(j=0;j<nfp;j++){ // loop over vertices
        u[0]=eloc[e][1][d*j];
        u[1]=eloc[e][1][d*j+1];
 
        xvertA[2*j] = 0;
        xvertA[2*j+1] = 0;
        for(k=0;k<nbasis;k++){
          bv[k] = basis[e][k](u);
          xvertA[2*j]   += bv[k]*xcutA[ixc+k];
          xvertA[2*j+1] += bv[k]*xcutA[ixc+k+3];
        }
      }
if(debug){
printf("\tcut vertices at: [%f %f; %f %f; %f %f]\n",xvertA[0],xvertA[1],xvertA[2],xvertA[3],xvertA[4],xvertA[5]);
}
      // Find the bounds of the overset boundary
      if(isBetween(xvertA[0],xvertA[1],xvertA[2],xvertA[3],pt[0],pt[1])){
	xseg[0] = xvertA[0];
	xseg[1] = xvertA[1];
	xseg[2] = xvertA[2];
	xseg[3] = xvertA[3];
      } else if(isBetween(xvertA[0],xvertA[1],xvertA[4],xvertA[5],pt[0],pt[1])){
	xseg[0] = xvertA[0];
	xseg[1] = xvertA[1];
	xseg[2] = xvertA[4];
	xseg[3] = xvertA[5];
      } else {
	xseg[0] = xvertA[2];
	xseg[1] = xvertA[3];
	xseg[2] = xvertA[4];
	xseg[3] = xvertA[5];
      }
if(debug){
printf("\toverset edge at: [%f %f; %f %f]\n",xseg[0],xseg[1],xseg[2],xseg[3]);
}

      // break up overset face by element 
      // fill OSFwgt and OSFshpL
      createOversetSegs(xseg,JinvA+iptrA[pc*eid+8],necutB,xcutB,cut2eB,iptrcB,OSFnseg+i,OSFwgt+iosf,OSFshpL+ishp,pccut,e,p,d,debug);
if(debug) printf("\tnseg = %i\n",OSFnseg[i]);

      // fill OSFshpR 
//      interpOversetCutNodes(xA+ix, xB, iptrB, pc, cutoversetA+ico, bfcutRA+ibf, JinvB, d, e, p, nelemB,debug);  // get cell id and shp value from mesh B

// need to empty bfcutL and R for overset faces
    }
  }

  // free pointers
  dgsand_free(xseg); 
//  dgsand_free(xgseg); 

  exit(1);       
}

