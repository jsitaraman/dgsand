void JacP1Tri( double *jac, double *xtmp, double* det){
//      N1 = 1-r-s, N2 = r, N3 = s
//          J(1,1) = x2-x1
//          J(1,2) = x3-x1
//          J(2,1) = y2-y1
//          J(2,2) = y3-y1
  jac[0] = xtmp[1]-xtmp[0];
  jac[1] = xtmp[2]-xtmp[0];
  jac[2] = xtmp[4]-xtmp[3];
  jac[3] = xtmp[5]-xtmp[3];
  *det = jac[0]*jac[3]-jac[1]*jac[2];
}

int FIND_NECUT(double x0, double *x,int* iptr, int d, int e, int p, int nelem, int pc, int imesh)
{
//Hack just to find necut. Only need this for debugging cases, 
//We just need necut so we can allocate all the other matrices
//
  // Define the cutting boundary as x = x0 (cutting away x < x0)

  int nfp=facePerElem[e];
  int i,f,j,k,l,m,n,tally[nfp],sum,ix,jp1,jp2,j0,fid;

  double jac[4],det=-1000; 
  double u[2],a,b,ycut1,ycut2;
  double xvert[6]; // x1 y1 x2 y2 x3 y3
  double xtmp[6];
  int vcut[2],vorig[2],neigh[3],forig[3];

  int nbasis=order2basis[e][p]; // first order bases are interpolative at nodes so it's ok to do this?
  double bv[nbasis*nfp];
  int necut = 0;
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
	if(imesh==0) if(xvert[2*j]>x0) tally[j]=1;
	if(imesh==1) if(xvert[2*j]<x0) tally[j]=1;
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

      //increment necut by 1 or 2 depending on how it was cut
      if(sum<nfp) necut+=sum; 
  }
printf("NECUT = %i\n",necut); 
  return necut; 
}

void CUT_CELLS(double x0, double *x, double* xcut, int* iptr, int* cut2e, int d, int e, int p, int nelem, int pc, int *cut2face, int* cut2neigh, int* elem2face, int* faces, int* iblank, int* cutoverset, int imesh, int ng, int* cellmerge)
// This routine cuts the cells according to some arbitrary vertical line. 
// This is for testing purposes and will eventually be replaced with 
// an actual cutting routine.
//
//OUTPUTS: 
//  xcut - global node coords of the cut cells
//  cut2e - map between cut cell id and original elem id
//  cut2neigh - map between cut face and R side element neighbor
//  cutoverset - tag list of cell faces that have overset boundaries
//  cellmerge - 0: cell uncut, 1: cell cut but not severly, 2: cell cut and needs merging
{

  // Define the cutting boundary as x = x0 (cutting away x < x0)

  int nfp=facePerElem[e];
  int i,f,j,k,l,m,n,tally[nfp],sum,ix,jp1,jp2,j0,fid;

  double jac[4],det=-1000; 
  double u[2],a,b,ycut1,ycut2;
  double xvert[6]; // x1 y1 x2 y2 x3 y3
  double xtmp[6];
  int vcut[3],vorig[3],neigh[3],forig[3];

  int nbasis=order2basis[e][p]; // first order bases are interpolative at nodes so it's ok to do this?
  double bv[nbasis*nfp];
  double area=0, area_cut=0;
//  int ng = ngGL[e][p]; 

  n = 0; 
  for(i=0;i<nelem;i++){
      iblank[i] = 0; 
      cellmerge[i] = 0;
      ix=iptr[pc*i+1];
      // get bases at rst = [0 0 ; 1 0; 0 1]
      for(j=0;j<nfp;j++){ // loop over vertices
        tally[j]=0;
        u[0]=eloc[e][1][d*j];
        u[1]=eloc[e][1][d*j+1];

        // get the vertex global coordinates of the original triangle
        // by doing x = sum(Ni * xi)
	xvert[2*j] = 0;
	xvert[2*j+1] = 0;
        for(k=0;k<nbasis;k++){ // accumulate bases and vertex coords
          bv[j*nbasis+k] = basis[e][k](u);
          xvert[2*j]   += bv[j*nbasis+k]*x[ix+k];
          xvert[2*j+1] += bv[j*nbasis+k]*x[ix+k+nbasis];
        }         

	//keep track of how many vertices are on cut side of cut boundary
	if(imesh==0) if(xvert[2*j]>x0) tally[j]=1;
	if(imesh==1) if(xvert[2*j]<x0) tally[j]=1;
      }

      // compute area of complete triangle
      // need to reorder xvert b/c I wrote this is in a stupid way
      // and changed convention
      xtmp[0] = xvert[0];
      xtmp[1] = xvert[2];
      xtmp[2] = xvert[4];
      xtmp[3] = xvert[1];
      xtmp[4] = xvert[3];
      xtmp[5] = xvert[5];
      JacP1Tri(jac,xtmp,&area); 
      area = 0.5*fabs(area); 
      area_cut = 0.0; 

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
      else if(sum==nfp){
        iblank[i] = 1;
        printf("\nelem %i is blanked\n",i+1);
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
          if(vcut[j] >= 0){
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
  	xtmp[3] = xvert[2*j0+1];
        xtmp[1] = x0;
	xtmp[4] = ycut1;
        xtmp[2] = x0;
	xtmp[5] = ycut2;

	// get the element neighbors on edges 0 and 2

        fid = abs(elem2face[i*3+j0])-1; //track which global face was cut
        forig[0] = fid; 
        neigh[0] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 

        forig[1] = -1; 
	neigh[1] = -2; // overset boundary

        fid = abs(elem2face[i*3+jp2])-1;
        forig[2] = fid; 
        neigh[2] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 

 	// initialize cutoverset index array
 	for(int a=0;a<3;a++) cutoverset[3*n+a] = -1; 

        // double check for positive jacobian
        // and save out info in correct order
        JacP1Tri(jac,xtmp,&det); 
//printf("%i cut area = %f\n",i,fabs(det)/2); 
        area_cut += fabs(det)/2;
        if(det>0){
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[1];
          xcut[n*3*2+2] = xtmp[2];

          xcut[n*3*2+3] = xtmp[3];
          xcut[n*3*2+4] = xtmp[4];
          xcut[n*3*2+5] = xtmp[5];

	  cut2face[n*3]   = forig[0];
	  cut2face[n*3+1] = forig[1];
	  cut2face[n*3+2] = forig[2];

	  cut2neigh[3*n] = neigh[0];
	  cut2neigh[3*n+1] = neigh[1];
	  cut2neigh[3*n+2] = neigh[2];

          cutoverset[3*n + 1] = 1; 
//	  for(int a=0;a<ng;a++) cutoverset[3*ng*n + 1*ng + a] = 1; 
        }
        else{
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[2];
          xcut[n*3*2+2] = xtmp[1];

          xcut[n*3*2+3] = xtmp[3];
          xcut[n*3*2+4] = xtmp[5];
          xcut[n*3*2+5] = xtmp[4];

	  cut2face[n*3]   = forig[0];
	  cut2face[n*3+1] = forig[2];
	  cut2face[n*3+2] = forig[1];

	  cut2neigh[3*n] = neigh[0];
	  cut2neigh[3*n+1] = neigh[2];
	  cut2neigh[3*n+2] = neigh[1];

          cutoverset[3*n + 2] = 1; 
//	  for(int a=0;a<ng;a++) cutoverset[3*ng*n + 2*ng + a] = 1; 
        }
        cut2e[n]=i; //store id of orig elem

	printf("\norig elem %i, cut cell %i, det %f:\n",i,n,det);
	if(fabs(det)<1e-13){
	printf("\nERROR: det %f %f \n jac = %f %f %f %f\nxtmp = %f %f %f %f %f %f\n",det,fabs(det),jac[0],jac[1],jac[2],jac[3],xtmp[0],xtmp[1],xtmp[2],xtmp[3],xtmp[4],xtmp[5]);
	for(f=0;f<nfp;f++){
	for(j=0;j<nbasis;j++){
	printf("f = %i,  bv[%i] = %f\n",f,j,bv[f*nbasis+j]);
	}
	}
	}

	printf("\tfull vtxcoords: (%f, %f), (%f, %f), (%f, %f)\n", xvert[0],xvert[1],xvert[2],xvert[3],xvert[4],xvert[5]);
	printf("\tcut vtx coords: (%f, %f), (%f, %f), (%f, %f)\n",xcut[n*3*2], xcut[n*3*2+3], xcut[n*3*2+1], xcut[n*3*2+4], xcut[n*3*2+2], xcut[n*3*2+5]);
	printf("\tcut elem neigh(%i, %i, %i) = %i %i %i\n",3*n,3*n+1,3*n+2,cut2neigh[3*n],cut2neigh[3*n+1],cut2neigh[3*n+2]);
	printf("\tcut2face = %i %i %i\n",cut2face[3*n],cut2face[3*n+1],cut2face[3*n+2]);
        n++; 
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
        xtmp[3] = xvert[2*j0+1];
        xtmp[1] = xvert[2*jp1];
        xtmp[4] = xvert[2*jp1+1];
        xtmp[2] = x0; 
        xtmp[5] = ycut1;

        fid = abs(elem2face[i*3+j0])-1;
        forig[0] = fid; 
        neigh[0] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 

        fid = abs(elem2face[i*3+jp1])-1;
        forig[1] = fid; 
        neigh[1] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 

	forig[2] = -1; 
	neigh[2] = i;  // this face's fluxes cancel out

 	// initialize cutoverset index array
 	for(int a=0;a<3;a++) cutoverset[3*n+a] = -1; 

        // double check for positive jacobian 
        JacP1Tri(jac,xtmp,&det); 
        area_cut += fabs(det)/2;
        if(det>0){
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[1];
          xcut[n*3*2+2] = xtmp[2];
          xcut[n*3*2+3] = xtmp[3];
          xcut[n*3*2+4] = xtmp[4];
          xcut[n*3*2+5] = xtmp[5];
          
	  cut2face[3*n] = forig[0];
	  cut2face[3*n+1] = forig[1];
	  cut2face[3*n+2] = forig[2];

	  cut2neigh[3*n] = neigh[0];
	  cut2neigh[3*n+1] = neigh[1];
	  cut2neigh[3*n+2] = neigh[2];        
        }
        else{
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[2];
          xcut[n*3*2+2] = xtmp[1];
          xcut[n*3*2+3] = xtmp[3];
          xcut[n*3*2+4] = xtmp[5];
          xcut[n*3*2+5] = xtmp[4];

	  cut2face[3*n] = forig[0];
	  cut2face[3*n+1] = forig[2];
	  cut2face[3*n+2] = forig[1];

	  cut2neigh[3*n] = neigh[0];
	  cut2neigh[3*n+1] = neigh[2];
	  cut2neigh[3*n+2] = neigh[1];
        }
 
        cut2e[n]=i; //store id of orig elem
	printf("\norig elem %i, cut cell %i, det %f:\n",i,n,det);
	if(fabs(det)<1e-13){
	printf("\nERROR: det %f %f \n jac = %f %f %f %f\nxtmp = %f %f %f %f %f %f\n",det,fabs(det),jac[0],jac[1],jac[2],jac[3],xtmp[0],xtmp[1],xtmp[2],xtmp[3],xtmp[4],xtmp[5]);
	for(f=0;f<nfp;f++){
	for(j=0;j<nbasis;j++){
	printf("f = %i,  bv[%i] = %f\n",f,j,bv[f*nbasis+j]);
	}
	}
	}
	printf("\tfull vtxcoords: (%f, %f), (%f, %f), (%f, %f)\n", xvert[0],xvert[1],xvert[2],xvert[3],xvert[4],xvert[5]);
	printf("\tcut vtx coords: (%f, %f), (%f, %f), (%f, %f)\n",xcut[n*3*2], xcut[n*3*2+3], xcut[n*3*2+1], xcut[n*3*2+4], xcut[n*3*2+2], xcut[n*3*2+5]);
	printf("\tcut elem neigh(%i, %i, %i) = %i %i %i\n",3*n,3*n+1,3*n+2,cut2neigh[3*n],cut2neigh[3*n+1],cut2neigh[3*n+2]);
	printf("\tcut2face = %i %i %i\n",cut2face[3*n],cut2face[3*n+1],cut2face[3*n+2]);
        n++; // end of first triangle

	//Second triangle (cut node 1, both intersect nodes)
        xtmp[0] = x0; 
        xtmp[3] = ycut1; 
        xtmp[1] = x0; 
        xtmp[4] = ycut2;
        xtmp[2] = xvert[2*j0]; 
        xtmp[5] = xvert[2*j0+1]; 

/*	neigh[0] = -2; // this is the overset boundary
        fid = abs(elem2face[i*3+jp2])-1;
        neigh[1] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2]; 
	neigh[2] = i; // this edge's fluxes cancel out 
*/

        fid = abs(elem2face[i*3+j0])-1;
        forig[0] = -i; // overset boundary
        neigh[0] = -2; // overset boundary 

        fid = abs(elem2face[i*3+jp2])-1;
        forig[1] = fid;
        neigh[1] = faces[6*fid+2] == i ? faces[6*fid+4] : faces[6*fid+2];

        fid = abs(elem2face[i*3+jp2])-1;
        forig[2] = -1; // this edge's fluxes cancel out
        neigh[2] = i; 
	
 	// initialize cutoverset index array
 	for(int a=0;a<3;a++) cutoverset[3*n+a] = -1; 

        // double check for positive jacobian 
        JacP1Tri(jac,xtmp,&det); 
        area_cut += fabs(det)/2;
        if(det>0){
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[1];
          xcut[n*3*2+2] = xtmp[2];
          xcut[n*3*2+3] = xtmp[3];
          xcut[n*3*2+4] = xtmp[4];
          xcut[n*3*2+5] = xtmp[5];

          cut2face[3*n] = forig[0];
          cut2face[3*n+1] = forig[1];
          cut2face[3*n+2] = forig[2];

	  cut2neigh[3*n] = neigh[0];
	  cut2neigh[3*n+1] = neigh[1];
	  cut2neigh[3*n+2] = neigh[2];
        }
        else{
          xcut[n*3*2]   = xtmp[0];
          xcut[n*3*2+1] = xtmp[2];
          xcut[n*3*2+2] = xtmp[1];
          xcut[n*3*2+3] = xtmp[3];
          xcut[n*3*2+4] = xtmp[5];
          xcut[n*3*2+5] = xtmp[4];

          cut2face[3*n] = forig[0];
          cut2face[3*n+1] = forig[2];
          cut2face[3*n+2] = forig[1];

	  cut2neigh[3*n] = neigh[0];
	  cut2neigh[3*n+1] = neigh[2];
	  cut2neigh[3*n+2] = neigh[1];
        } 
	cutoverset[3*n + 0] = 1; 

        cut2e[n]=i; //store id of orig elem
	printf("\norig elem %i, cut cell %i, det %f:\n",i,n,det);
	if(fabs(det)<1e-13){
	printf("\nERROR: det %f %f \n jac = %f %f %f %f\nxtmp = %f %f %f %f %f %f\n",det,fabs(det),jac[0],jac[1],jac[2],jac[3],xtmp[0],xtmp[1],xtmp[2],xtmp[3],xtmp[4],xtmp[5]);
	for(f=0;f<nfp;f++){
	for(j=0;j<nbasis;j++){
	printf("f = %i,  bv[%i] = %f\n",f,j,bv[f*nbasis+j]);
	}
	}
	}
	printf("\tfull vtxcoords: (%f, %f), (%f, %f), (%f, %f)\n", xvert[0],xvert[1],xvert[2],xvert[3],xvert[4],xvert[5]);
	printf("\tcut vtx coords: (%f, %f), (%f, %f), (%f, %f)\n",xcut[n*3*2], xcut[n*3*2+3], xcut[n*3*2+1], xcut[n*3*2+4], xcut[n*3*2+2], xcut[n*3*2+5]);
	printf("\tcut elem neigh(%i, %i, %i) = %i %i %i\n",3*n,3*n+1,3*n+2,cut2neigh[3*n],cut2neigh[3*n+1],cut2neigh[3*n+2]);
	printf("\tcut2face = %i %i %i\n",cut2face[3*n],cut2face[3*n+1],cut2face[3*n+2]);
        n++; // end of second triangle
      } // if cut cell
      // Check if cell needs merging
      if(area_cut<1e-14){
	cellmerge[i] = 0;
      }
      else if(area_cut/area < 0.5){
        cellmerge[i] = 1; // cut but doesn't require cell merging
      }
      else{
        cellmerge[i] = 2; // requires cell merging
      }
      printf("\n\telem %i, area_cut = %f, area = %f, cellmerge = %i\n",i,area_cut,area,cellmerge[i]);
  } // nelem loop

  // Debug
  for(int aa = 0;aa<n;aa++)
    for(int f = 0;f<nfp;f++){
      printf("\norig %i, el %i, side %i = %i\n",cut2e[aa],aa,f);
      printf("\tcutoverset[%i] = %i\n",aa*nfp+f,cutoverset[aa*nfp+f]);
    }
}

