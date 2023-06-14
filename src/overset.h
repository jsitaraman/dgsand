#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void setOversetFluxes(double *OSFflux, int OSFnseg, int* OSFeID, double* OSFshpL, double* OSFshpR, double* qA, double* qB, int* iptrB, int* elemParentB, int d, int e, int p, int pc, int pde, int debug)
// computes R q values at cut cell interface
// using q from other mesh
{
  int nfp = facePerElem[e];
  int nbasis = order2basis[e][p];
  int ngauss = ngGL[e][p]; 
  int nfields=get_nfields[pde](d);
  int eid, pid, iq;   

  if(debug) printf("in setOversetFluxes\n"); 

  // loop over segments
  int floc = 0; 
  int bloc = 0; 
  for(int i=0;i<OSFnseg;i++){
    for(int j=0;j<ngauss;j++){
      eid = OSFeID[i*ngauss+j];
      pid = elemParentB[eid];
      iq = iptrB[pc*pid]; 

      if(debug) printf("seg %i, gauss %i\n",i,j); 
      if(debug) for(int b=0;b<nbasis;b++) printf("\t\tdebug: OSFshpL[%i] = %f, OSFshpR[%i] = %f\n",
                                                  bloc+b,OSFshpR[bloc+b],bloc+b,OSFshpL[bloc+b]);

      for(int f=0;f<nfields;f++){
        // L state
	OSFflux[floc + f] = OSFshpL[bloc]*qA[f*nbasis];
	for(int b=1;b<nbasis;b++) 
	  OSFflux[floc + f] += OSFshpL[bloc+b]*qA[f*nbasis+b];

        // R state 
	OSFflux[floc + nfields + f] = OSFshpR[bloc]*qB[iq+f*nbasis];
	for(int b=1;b<nbasis;b++) 
	  OSFflux[floc + nfields + f] += OSFshpR[bloc+b]*qB[iq+f*nbasis+b];	  
      }// loop over fields

      if(debug){
        for(int f=0;f<nfields;f++)
        for(int b=0;b<nbasis;b++) 
          printf("qA[%i] = %f, qB[%i] = %f\n",f*nbasis+b,qA[f*nbasis+b],
                 iq+f*nbasis+b,qB[iq+f*nbasis+b]);
        for(int f=0;f<nfields;f++)
          printf("L state = %f, R state = %f\n",OSFflux[floc + f],OSFflux[floc + nfields + f]);
      }
      for(int aa = 0; aa<3*nfields; aa++){
        if(isnan(OSFflux[floc+aa])){
          printf("\nERROR: nan in OSFflx:\n");
            for(int bb = 0; bb<3*nfields; bb++) printf("\tOSFflx[%i] = %.16e\n",bb,OSFflux[floc+bb]);
            exit(1);
        }
      } 

      floc = floc + 3*nfields;
      bloc = bloc + nbasis; 
    } // loop over gauss pts on each segment
  } // loop over segments
}

void EXCHANGE_OVERSET(double* OSFflux, double* OSFshpL, double* OSFshpR, int* OSFnseg, int* OSFeID, 
                      int* cut2eA, double* qA, double* qB, int* cutoverset, 
                      int* iptrcA, int* iptrA, int* iptrB, int* elemParentA, int* elemParentB, 
                      int necutA, int pccut, int d, int e, int p, int pc, int pde, int imesh)
{
  int iel,ishp, ip, ix, iq, ibf, ic2n, iflx, ico;
  int nfp = facePerElem[e];
  int eid, pid, flag, debug, m;
  int nbasis=order_to_basis(e,p);        // basis for solution
  // Loop over cut cells in mesh A
  for(int i = 0; i<necutA; i++){
    flag = 0; 
    ip=pccut*i;
    ico  = iptrcA[ip+12]; 
    iel  = iptrcA[ip+14]; 
    ishp = iptrcA[ip+15]; 
    iflx = iptrcA[ip+16];

    eid = cut2eA[i];
    pid = elemParentA[eid];
    iq = iptrA[pid*pc]; // using parent's q values

    if(i==2 && eid == 1){
      debug = 1;
    } else{
      debug=0;
    }
    if(debug) printf("In Exchange\n");

    // does this cut cell have an overset boundary?
    if(cutoverset[ico] + cutoverset[ico+1] + cutoverset[ico+2]>-nfp){
      // interpolate q fluxes from mesh B
      setOversetFluxes(OSFflux+iflx, OSFnseg[i],OSFeID+iel, OSFshpL+ishp, OSFshpR+ishp, 
                       qA+iq, qB, iptrB, elemParentB, d, e, p, pc, pde, debug);
    } // if cut overset
  } // cut cells A loop
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

int isBetween(double xA, double yA, double xB, double yB, double xC, double yC)
// checks to see if (xC,yC) is on the straight line between (xA,yA) and (xB,yB)
{
  double d1 = sqrt((xA-xC)*(xA-xC) + (yA-yC)*(yA-yC));
  double d2 = sqrt((xB-xC)*(xB-xC) + (yB-yC)*(yB-yC));
  double d3 = sqrt((xA-xB)*(xA-xB) + (yA-yB)*(yA-yB));

  double tol = 1e-10; 
//printf("isbetween: d = %f %f %f, abs = %f\n",d1,d2,d3,fabs(d3-d1-d2));
/*  if(fabs(d3-d1-d2)<tol && 
     (fabs(xC-xA)>tol && fabs(yC-yA)>tol) &&
     (fabs(xC-xB)>tol && fabs(yC-yB)>tol)){
printf("isBetween conditions: %.13e %i %i %i\n",fabs(d3-d1-d2),!(fabs(xC-xA)>tol && fabs(yC-yA)>tol),!(fabs(xC-xB)>tol && fabs(yC-yB)>tol),(fabs(d3-d1-d2)<tol &&      !(fabs(xC-xA)>tol && fabs(yC-yA)>tol) &&     !(fabs(xC-xB)>tol && fabs(yC-yB)>tol)));
*/

  if(fabs(d3-d1-d2)<tol){
    return 1;
  }
  else{
    return 0; 
  }

}

int isSame(double A, double B){
  double tol = 1e-12; 
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

  // Assuming first 2 entries into in array are x and y of one of the end pts
  double ds[n];
  out[0] = in[0];
  out[1] = in[1];

  // calculate distance s from first node
  ds[0] = 1e15;
  for(int i=1;i<n;i++)
    ds[i]= pow(in[2*i]-in[0],2) + pow(in[2*i+1]-in[1],2); 

  // loop over out
  double minS;
  int ind;
  for(int i = 1;i<n;i++){ // loop over out
    minS = 1e15;
    for(int j = 1;j<n;j++){ // loop over in
      if(ds[j]<minS){
        minS = ds[j];
        ind = j; 
      }
    } // loop over in
   
    ds[ind] = 1e15;
    out[2*i] = in[2*ind];
    out[2*i+1] = in[2*ind+1];
  }
}

void createOversetGauss(double* xA, double* xB, double* xseg, double* JinvA, double* JinvB, 
                        int necutB, double *xcutA, double* xcutB, int* cut2eB, 
                        int* iptrB, int* iptrcB, int* elemParentB,
			int* OSFeID, int* OSFnseg, 
                        double* OSFxn, double* OSFshpL, double* OSFshpR, 
                        int pc, int pccut,int e, int p, int d, int nelemB, 
                        int debug,int iface, int ixn, int ismerge)
// Break up overset face into segments, one for each overlapping element, and distribute quad pts
{
  int ij,ip, eid, ixc,ixB; 
  int nfp = facePerElem[e];
  int b,i,j,k,m,n,f,uniq;
  int nbasis=order_to_basis(e,p);        // basis for solution
  int nbasisx=order_to_basis(e,1);        // basis for solution
  double bd[nbasis][d];  // basis derivative  
  double mat[d][d], Ja[3];
  double Jb[3]={0,0,1};
  int maxseg = 20; 
  double u[2], bv[nbasis],pt[2],pt2[2],xvert[6];
  int ngGL=get_ngGL(e,p);
  double xcpy[maxseg*ngGL*d]; // tmp copy of xseg
  double x0[2],xN[2],xloc,yloc; 
  
  // get x0 for the full element
  x0[0] = 0.0;
  x0[1] = 0.0;
  u[0] = 0.0;
  u[1] = 0.0;
  for(b=0;b<nbasis;b++){
    x0[0] +=  basis[e][b](u)*xA[b]; 
    x0[1] +=  basis[e][b](u)*xA[b+nbasis];
  } 

  // copy overset end pt coords, rest should be blank anyways
  for(i=0;i<2*d;i++) xcpy[i]=xseg[i]; 

  // Figure out how many segments on face
  // loop over cut cells on mesh B
  OSFnseg[0] = 1; // all have at least 1 seg
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
      for(k=0;k<nbasisx;k++){
        bv[k] = basis[e][k](u);  
        pt[0] += bv[k]*xcutB[ixc+k];
        pt[1] += bv[k]*xcutB[ixc+k+3];
      }

      // if vert is on cut overset, then store mesh B original cell ID and vertex coordinate
      if(isBetween(xcpy[0],xcpy[1],xcpy[2],xcpy[3],pt[0],pt[1])){
        if(debug) printf("Found new pt at (%f %f), curr nseg = %i\n",pt[0],pt[1],OSFnseg[0]);

        // store vertex coordinate if unique
        uniq = 1;
        for(k=0;k<(OSFnseg[0]+1);k++){          
          if(debug) 
            printf("\t check unique. pt = (%f %f), ref = (%f %f)\n",pt[0],pt[1],xcpy[2*k],xcpy[2*k+1]);
          if(isSame(xcpy[2*k],pt[0]) && isSame(xcpy[2*k+1],pt[1])) uniq = 0;
        }
	
        if(uniq){
          OSFnseg[0]++;                   	  	  
          xcpy[2*(OSFnseg[0])]   = pt[0];
          xcpy[2*(OSFnseg[0])+1] = pt[1];
	  if(debug) printf("\tPoint is unique, now nseg = %i\n",OSFnseg[0]);
        }
      } // found vertex in between
    } // loop over edges B
  } // loop over cut cell B

  // debug
  if(debug){
    printf("\n\t%i segments. Unordered list: \n",OSFnseg[0]);
    for(i=0;i<OSFnseg[0]+1;i++)
      printf("\t\t%f %f\n",xcpy[2*i],xcpy[2*i+1]);
  }
  
  // Reorder list of segment end pts if necessary
  if(OSFnseg[0]>1){
    orderCoords(xcpy,xseg,OSFnseg[0]+1);
  } else{
    xseg=xcpy; 
    OSFnseg[0] = 1;
  }

  // debug
  if(debug){
  printf("\n\tOrdered list\n"); 
  for(i=0;i<OSFnseg[0]+1;i++)
    printf("\t\t%f %f\n",xseg[2*i],xseg[2*i+1]);
  }

  // Distribute gauss pts between each of the segments
  double nx, ny, s, L;
  int aa,bb,cc;
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

      xloc = pt[0] + nx*s;		// gauss pt xy coord
      yloc = pt[1] + ny*s;

      //=====================
      // Fill in OSFshpL and OSFxn
      //=====================
      // convert physical coord into parent element rst coord
      // Note: reusing pt2 here
      // rst = [JinvA][xy]
      pt2[0] = xloc-x0[0];
      pt2[1] = yloc-x0[1];
      // Using JinvA from parent element of msh A
      axb(JinvA,pt2,u,d);
      printf("===============\n");
      for(k=0;k<nbasis;k++){
	OSFshpL[m+k] = basis[e][k](u);
        if(ismerge==0 && OSFshpL[m+k]<0 && fabs(OSFshpL[m+k])>1e-14){
	  printf("ERROR: OSFshpL less than 0. Exiting\n");
	  printf("\tOSFshpL = %.16e\n",OSFshpL[m+k]);
	  printf("\txloc = %f %f\n",xloc,yloc);
	  printf("\tx0 = %f %f\n",x0[0],x0[1]);
	  printf("\tpt2 = %f %f\n",pt2[0],pt2[1]);
	  printf("\tJinvA = %f %f; %f %f\n",JinvA[0],JinvA[1],JinvA[2],JinvA[3]);
	  printf("\tu = %.16e %.16e\n",u[0],u[1]);
	  exit(1);
	}
        if(debug) printf("\t\t\tOSFshpL %i = %f\n",k,OSFshpL[m+k]);
      }

      for(b=0;b<nbasis;b++)
        for(k=0;k<d;k++)
          bd[b][k]=basis_d[e][b*d+k](u);

      // compute wall normals here
      // build jacobian
      for(aa=0;aa<d;aa++){
        Ja[aa]=0;
        for(bb=0;bb<d;bb++){
          mat[aa][bb]=xcutA[aa*3]*bd[0][bb];
          for(cc=1;cc<3;cc++)
            mat[aa][bb]+=xcutA[aa*3+cc]*bd[cc][bb];
        } // loop over d
        for(bb=0;bb<d;bb++)
          Ja[aa]+=(mat[aa][bb]*face2elem[e][d*iface+bb]); // face2elem goes from edge coord to rst elem coord
      } // aa loop over d
if(debug) printf("pre Ja = %f %f %f %f\n",Ja[0],Ja[1],Ja[2],Ja[3]);

      // scale Ja to be segment length
      double norm = 0.0;
      for(aa=0;aa<3;aa++)
	norm+=Ja[aa]*Ja[aa];
      norm=sqrt(norm);
      for(aa=0;aa<3;aa++)
	Ja[aa] = Ja[aa]*L/norm;

      // do Weighted normal = Ja x zhat
      // This gives me the normal vector
      // pointing outward of the cut triangle
      cross(&(OSFxn[i*d*ngGL+j*d]),Ja,Jb,d);

      if(isnan(OSFxn[d*i*ngGL + d*j]) || isnan(OSFxn[d*i*ngGL + d*j +1] ||
         fabs(OSFxn[d*i*ngGL + d*j]) + fabs(OSFxn[d*i*ngGL + d*j+1]) < 1e-12)){
	printf("\nERROR: normal %i is nan: %f %f\n",ixn+d*i*ngGL+d*j,OSFxn[d*i*ngGL + d*j], OSFxn[d*i*ngGL + d*j +1]);
	printf("Ja = %f %f %f\n",Ja[0],Ja[1],Ja[2]);
	printf("mat = %f %f %f %f\n",mat[0][0],mat[0][1],mat[1][0],mat[1][1]);
	printf("face2elem %f %f\n",face2elem[e][d*iface],face2elem[e][d*iface+1]);
	exit(1);
      }

      // Normal vector needs to be flipped 
      // to be outward from remaining triangle
      for(aa=0;aa<d;aa++) OSFxn[d*i*ngGL+d*j+aa] = -OSFxn[d*i*ngGL+d*j+aa];
      if(debug){
	printf("\t\tgauss pt %i at [%f %f] normal %i %f %f\n",i*ngGL+j,
	       xloc,yloc,ixn+d*i*ngGL+d*j,OSFxn[d*i*ngGL + d*j], OSFxn[d*i*ngGL + d*j +1]);
      }

      //=====================
      // Fill in OSFshpR
      //=====================

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
  	  OSFeID[i*ngGL+j] = elemParentB[n]; 	// store mesh B parent elemID for this overset gauss pt 
          break; // exit loop over mesh B
	}
      }

      // compute the rst coordinate
      pt2[0] = xloc-xvert[0];
      pt2[1] = yloc-xvert[1];
      ij = iptrB[pc*OSFeID[i*ngGL+j]+4]; 
      axb(JinvB+ij,pt2,u,2);

      // if finished loop and still not found, something is wrong
      if(inside==0 || (n==elemParentB[n] && (u[0]<0 || u[0]>1 || u[1]<0 || u[1]>1))){
        printf("ERROR! Can't find mesh B element for mesh A point (%f, %f)!\n",xloc,yloc); 
        printf("\tu = %f %f, inside = %i\n",u[0],u[1],inside);
        if(inside) printf("\tTri coords: (%f, %f), (%f, %f), (%f, %f)\n",
           xvert[0],xvert[1],xvert[2],xvert[3],xvert[4],xvert[5]);
        exit(0); 
      }
      else{
        if(debug){
          printf("\t\t\tfound mesh A pt (%f, %f) in mesh B elem %i (parent is %i): (%f, %f), (%f, %f), (%f, %f)\n",xloc,yloc,n,OSFeID[i*ngGL+j],
                 xvert[0],xvert[1],xvert[2],xvert[3],xvert[4],xvert[5]);
          printf("\t\t\trst = %f %f\n",u[0],u[1]);
        }
      }
      
      // get mesh B shape function values at quad pt and store in bfcutR
      for(int b=0; b<nbasis;b++){
        OSFshpR[m+b] = basis[e][b](u);
	if(debug)  printf("\t\t\tOSFshpR %i = %f\n",b,OSFshpR[m+b]);
      }

      m=m+nbasis; // counter for bases
    } // loop over gauss pts
  }  // loop over overset segments
}

void SETUP_OVERSET(int* cut2e, int* cut2eB, int* cutoversetA, 
		   int* iptrA, int* iptrB, int* iptrcA, int* iptrcB, 
		   double* xA, double* xB, double* xcutA, double* xcutB, 
		   double* bfcutLA, double* bfcutRA, double* JinvA, double* JinvB,
                   int* elemParentA, int* elemParentB, 
		   int* OSFnseg, int* OSFeID, double* OSFxn, double* OSFshpL, double* OSFshpR, 
		   int d, int e, int p, int pc, int pccut, 
		   int necutA, int necutB, int nelemB)
// For overset boundaries in mesh A, find corresponding elements and shape function values on mesh B. Build multiple segments along overset boundary to handle discontinuous overset flux across multiple elements
{
  int ip, ix, ixc, iq, ibf, ibfd, ifw, ic2n,iflx, cibf,flag, cid, ixn, iosf, ishp, ismerge;
  int nfp = facePerElem[e];
  int eid, pid, i,j,k,m, count,jp1,debug,iface;
  int nbasis=order_to_basis(e,p);        // basis for solution
  int nbasisx=order_to_basis(e,1);        // basis for solution
  int ngGL=get_ngGL(e,p);
  double pt[2],u[2],bv;
  int maxseg=20; 
  double *xseg;
  xseg=dgsand_alloc(double,((maxseg+1)*d)); // physical coordinates of segment end pts

  // Loop over cut cells in mesh A
  for(int i = 0; i<necutA; i++){
    ip=pccut*i;
    eid = cut2e[i]; 
    pid = elemParentA[eid];
    ix  =iptrA[eid*pc+1];
    ixc = iptrcA[ip+1]; 
    ibf =iptrcA[ip+6];
    ic2n=iptrcA[ip+12];
    ixn=iptrcA[ip+13];
    iosf=iptrcA[ip+14];
    ishp=iptrcA[ip+15];

    debug = 1;
    if(debug) printf("DEBUG cut elem %i, ixn = %i\n",i,ixn);

    // only handle elements with overset boundaries
    if(cutoversetA[ic2n]+cutoversetA[ic2n+1]+cutoversetA[ic2n+2]>-nfp){
      // get end pts of overset edge
      for(int j=0;j<nfp;j++){
        if(cutoversetA[ic2n+j]==1){
          jp1 = j+1 % 3;

	  // starting pt of current face
          u[0]=eloc[e][1][d*j];
          u[1]=eloc[e][1][d*j+1]; 
          xseg[0] = 0;
          xseg[1] = 0;
          for(k=0;k<nbasisx;k++){
            bv = basis[e][k](u);
            xseg[0] += bv*xcutA[ixc+k];  
            xseg[1] += bv*xcutA[ixc+k+3];
          }

	  // starting pt of next face
          u[0]=eloc[e][1][d*jp1];
          u[1]=eloc[e][1][d*jp1+1]; 
          xseg[2] = 0;
          xseg[3] = 0;
          for(k=0;k<nbasisx;k++){
            bv = basis[e][k](u);
            xseg[2] += bv*xcutA[ixc+k];
            xseg[3] += bv*xcutA[ixc+k+3];
          }
          iface = j;
        
	  if(isnan(xseg[0]) || isnan(xseg[1]) ||isnan(xseg[2]) ||isnan(xseg[3])){
	    printf("ERROR: nan in xseg\n");
	    printf("\t xseg = %.16e %.16e %.16e %.16e\n",xseg[0],xseg[1],xseg[2],xseg[3]);
	    for(k=0;k<nbasis;k++) 
	      printf("\t\t k = %i, xcutA = %.16e %.16e\n",k,xcutA[ixc+k],xcutA[ixc+k+3]);
	    exit(1); 
	  }

	  if(debug) 
	    printf("cut elem %i, face %i, xseg = [%f %f; %f %f]\n",i,j,xseg[0],xseg[1],xseg[2],xseg[3]);
        } // if overset
      } // loop faces

      // Distribute Gauss pts and fill OSFxn and OSFshpL and OSFshpR
      // Note: JinvA corresponds to parent elem, NOT current elem
      ismerge = (eid!=pid);
      createOversetGauss(xA+ix,xB,xseg,JinvA+iptrA[pc*pid+4],JinvB,necutB,
		         xcutA+ixc,xcutB,cut2eB,iptrB,iptrcB, elemParentB,
			 OSFeID+iosf,OSFnseg+i,OSFxn+ixn,OSFshpL+ishp,OSFshpR+ishp,
			 pc,pccut,e,p,d,nelemB,debug,iface,ixn,ismerge);
    }
  }

  // free pointers
  dgsand_free(xseg); 
}

