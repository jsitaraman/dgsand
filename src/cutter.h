void JacP1Tri( double *jac, double *xtmp, double det){
  jac[0] = xtmp[1]-xtmp[0];
  jac[1] = xtmp[2]-xtmp[0];
  jac[2] = xtmp[4]-xtmp[3];
  jac[3] = xtmp[5]-xtmp[3];
  det = jac[0]*jac[3]-jac[1]*jac[2];
}

void CUT_CELLS(double *x, double* xcut, int* iptr, int* iptrc, int* cut2e, int *necut, int d, int e, int nelem, int ncfaces, int pc, int* cut2neigh, int* elem2face, int* faces)
// This routine cuts the cells according to some arbitrary vertical line. 
// This is for testing purposes and will eventually be replaced with 
// an actual cutting routine.
//
//OUTPUTS: 
//  xcut - global node coords of the cut cells
//  cut2e - map between cut cell id and original elem id
//  necut - number of cut regions
//  ncfaces - number of faces of cut cells
//  /cut2neigh - map between cut face and R side element neighbor
{

// XXX I don't yet understand what ncfaces is supposed to mean.
// Need to verify this!

  // Define the cutting boundary as x = x0 (cutting away x < x0)
  double x0 = -0.5; // x coord where cut done

  int nfp=facePerElem[e];
  int i,j,k,l,m,n,tally[nfp],sum,ix,jp1,jp2,j0,neigh[3],fid;

  double jac[4],det; 
  double u[2],a,b,ycut1,ycut2;
  double xvert[6]; // x1 y1 x2 y2 x3 y3
  double xtmp[6];
  int vcut[2],vorig[2];

  int nbasis=order2basis[e][1]; // first order bases are interpolative at nodes so it's ok to do this?
  double bv[nbasis*nfp];
  (*necut) = 0;
  n = 0; 
  for(i=0;i<nelem;i++){
      ix=iptr[pc*i+1];
      // get bases at rst = [0 0 ; 1 0; 0 1]
      for(j=0;j<nfp;j++){ // loop over vertices
        tally[j]=0;
        u[0]=eloc[e][1][d*j];
        u[1]=eloc[e][1][d*j+1];

        // get the vertex global coordinates of the original trying
        // by doing x = sum(Ni * xi)
	xvert[2*j] = 0;
	xvert[2*j+1] = 0;
        for(k=0;k<nbasis;k++){ // accumulate bases and vertex coords
          bv[j] = basis[e][k](u);
          xvert[2*j]   += bv[j]*x[ix+k];
          xvert[2*j+1] += bv[j]*x[ix+k+nbasis];
        }         

	//keep track of how many vertices are on cut side of cut boundary
        if(xvert[2*j]<x0) tally[j]=1;
      }

      // check if element has vertices on both sides of x=x0 then 
      // track which vertices will be cut and which kept
      sum=0;
      for(j=0;j<nfp;j++) {
        sum+=tally[j];
        vcut[j]=-1;
        vorig[j]=-1;
      }
      if(sum>0 && sum<nfp){
        for(j=0;j<nfp;j++){ 
          if(tally[j]==1){
            vcut[j] = j;
          } else{
            vorig[j] = j;
          }
        }
      }

      // based on vorig and vcut, I can figure out which faces were cut 
      // and then fill out cut2neigh

      // need to store xcut (cut cell corner locations)
      // Currently only works for p=1 since that's all we need
      //
      // tri cut region (1 original node, 2 intersect pts)
      if(sum==1){ 
	for(j=0;j<nfp;j++){
	  neigh[j]=-1; 
          if(vcut[j] > 0){
            j0 = j;	// this is the cut node index
          }
        }
	jp1=(j0+1)%nfp;
	jp2=(j0+2)%nfp;

	// find the 2 intersection points
        a = (xvert[2*j0+1] - xvert[2*jp1+1])/(xvert[2*j0] - xvert[2*jp1]+1e-15);
        b = xvert[2*j0+1]-a*xvert[2*j0]; 
        ycut1=a*x0+b; 

        a = (xvert[2*jp2+1] - xvert[2*j0+1])/(xvert[2*jp2] - xvert[2*j0]+1e-15);
        b = xvert[2*j0+1]-a*xvert[2*j0]; 
        ycut2=a*x0+b; 

  	xtmp[0] = xvert[2*j0];
  	xtmp[1] = xvert[2*j0+1];
        xtmp[2] = x0;
	xtmp[3] = ycut1;
        xtmp[4] = x0;
	xtmp[5] = ycut2;

	// get the element neighbors on edges 0 and 2

        fid = abs(elem2face[i*3+j0])-1; //track which global face was cut
        neigh[0] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 
	neigh[1] = i;
        fid = abs(elem2face[i*3+jp2])-1;
        neigh[2] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 

        // double check for positive jacobian
        // and save out info in correct order
        JacP1Tri(jac,xtmp,det); 
        if(det>0){
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[2];
          xcut[n*3*2+2] = xtmp[4];

          xcut[n*3*2+3] = xtmp[1];
          xcut[n*3*2+4] = xtmp[3];
          xcut[n*3*2+5] = xtmp[5];
         
	  cut2neigh[n] = neigh[0];
	  cut2neigh[n+1] = neigh[1];
	  cut2neigh[n+2] = neigh[2];
        }
        else{
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[4];
          xcut[n*3*2+2] = xtmp[2];

          xcut[n*3*2+3] = xtmp[1];
          xcut[n*3*2+4] = xtmp[5];
          xcut[n*3*2+5] = xtmp[3];

	  cut2neigh[n] = neigh[0];
	  cut2neigh[n+1] = neigh[2];
	  cut2neigh[n+2] = neigh[1];
        }
        cut2e[n]=i; //store id of orig elem
        n++; 
        (*necut)++;
        ncfaces = ncfaces+3; 
      }

      // quad cut region (2 cut nodes, 2 intersect pts)
      // need to form 2 triangles here
      else if(sum==2){ 
	j0 = -1000; 
        for(j=0;j<nfp;j++){
          jp1=(j+1)%nfp;
	  if(vcut[j] >= 0 && vcut[jp1] >= 0){
            j0 = j; 
	  }
        }
	jp1 = (j0+1)%nfp; 
	jp2 = (j0+2)%nfp; 

	// get both intersect points	   
        a = (xvert[2*jp2+1] - xvert[2*jp1+1])/(xvert[2*jp2] - xvert[2*jp1]+1e-15);
        b = xvert[2*jp2+1]-a*xvert[2*jp2]; 
        ycut1=a*x0+b; 

        a = (xvert[2*j0+1] - xvert[2*jp2+1])/(xvert[2*j0] - xvert[2*jp2]+1e-15);
        b = xvert[2*j0+1]-a*xvert[2*j0]; 
        ycut2=a*x0+b; 

	// form first triangle (cut node 1 -> cut node 2 -> int 1)
        xtmp[0] = xvert[2*j0];
        xtmp[1] = xvert[2*j0+1];
        xtmp[2] = xvert[2*jp1];
        xtmp[3] = xvert[2*jp1+1];
        xtmp[4] = x0; 
        xtmp[5] = ycut1;

        fid = abs(elem2face[i*3+j0])-1;
        neigh[0] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 
        fid = abs(elem2face[i*3+jp1])-1;
        neigh[1] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 
	neigh[2] = i; 

        // double check for positive jacobian 
        JacP1Tri(jac,xtmp,det); 
        if(det>0){
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[2];
          xcut[n*3*2+2] = xtmp[4];
          xcut[n*3*2+3] = xtmp[1];
          xcut[n*3*2+4] = xtmp[3];
          xcut[n*3*2+5] = xtmp[5];
          
	  cut2neigh[n] = neigh[0];
	  cut2neigh[n+1] = neigh[1];
	  cut2neigh[n+2] = neigh[2];
        }
        else{
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[4];
          xcut[n*3*2+2] = xtmp[2];
          xcut[n*3*2+3] = xtmp[1];
          xcut[n*3*2+4] = xtmp[5];
          xcut[n*3*2+5] = xtmp[3];

	  cut2neigh[n] = neigh[0];
	  cut2neigh[n+1] = neigh[2];
	  cut2neigh[n+2] = neigh[1];
        }
 
        cut2e[n]=i; //store id of orig elem
        n++; 
        (*necut)++;
        ncfaces = ncfaces+3; 

	//Second triangle (cut node 1, both intersect nodes)
        xtmp[0] = x0; 
        xtmp[1] = ycut1; 
        xtmp[2] = x0; 
        xtmp[3] = ycut2;
        xtmp[4] = xvert[2*j0]; 
        xtmp[5] = xvert[2*j0+1]; 

	neigh[0] = i; // this is the overset boundary
        fid = abs(elem2face[i*3+jp2])-1;
        neigh[1] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 
	neigh[2] = i; // this edge's neighbor is itself
	
        // double check for positive jacobian 
        JacP1Tri(jac,xtmp,det); 
        if(det>0){
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[2];
          xcut[n*3*2+2] = xtmp[4];
          xcut[n*3*2+3] = xtmp[1];
          xcut[n*3*2+4] = xtmp[3];
          xcut[n*3*2+5] = xtmp[5];

	  cut2neigh[n] = neigh[0];
	  cut2neigh[n+1] = neigh[1];
	  cut2neigh[n+2] = neigh[2];
        }
        else{
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[4];
          xcut[n*3*2+2] = xtmp[2];
          xcut[n*3*2+3] = xtmp[1];
          xcut[n*3*2+4] = xtmp[5];
          xcut[n*3*2+5] = xtmp[3];

	  cut2neigh[n] = neigh[0];
	  cut2neigh[n+1] = neigh[2];
	  cut2neigh[n+2] = neigh[1];
        }
 
        cut2e[n]=i; //store id of orig elem
        n++; 
        (*necut)++;
        ncfaces = ncfaces+3; 
      }
  }

//store the xcut variables correctly using iptrc
  for(i=0;i<(*necut);i++){
    ix = iptrc[i*pc];
  }

}
