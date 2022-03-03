void checkJac( double *xtmp,double det){
  double jac[4]; 
  jac[0] = xtmp[0]-xtmp[2];
  jac[1] = xtmp[1]-xtmp[2];
  jac[2] = xtmp[3]-xtmp[5];
  jac[3] = xtmp[4]-xtmp[5];
  det = jac[0]*jac[3]-jac[1]*jac[2];
}

void CUT_CELLS(double *x, double* xcut, int* cut2e, int necut, int d, int e, int nelem, int pc)
// This routine cuts the cells according to some arbitrary vertical line. 
// This is for testing purposes and will eventually be replaced with 
// an actual cutting routine.
//
//OUTPUTS: 
//  xcut - global node coords of the cut cells
//  cut2e - map between cut cell id and original elem id
//  necut - number of cut regions
{

  // Define the cutting boundary as x = x0 (cutting away x < x0)
  double x0 = -0.5; // x coord where cut done

  int nfp=facePerElem[e];
  int i,j,k,l,m,n,tally[nfp],sum,ix;

  double jac[4],det; 
  double u[2],a,b,ycut1,ycut2;
  double xvert[6]; // x1 x2 x3 y1 y2 y3 
  double xtmp[6];
  double rst[6]={0,1,0,0,0,1}; // r1 r2 r3 s1 s2 s3

  int vcut[2],vorig[2];

  int nbasis=order2basis[e][1];
  double bv[nbasis*nfp];
  necut = 0;
  n = 0; 
  for(i=0;i<nelem;i++){
      ix=pc*i;

      // get bases at rst = [0 0 ; 1 0; 0 1]
      for(j=0;j<nfp;j++){ // loop over vertices
        tally[j]=0;
        u[0] = rst[j];
        u[1] = rst[j+3];
	xvert[2*j] = 0;
	xvert[2*j+3] = 0;

        // get the vertex global coordinates of the original trying
        // by doing x = sum(Ni * xi)
        for(k=0;k<nbasis;k++){ // accumulate bases and vertex coords
          bv[j] = basis[e][k](u);
          xvert[2*j] += bv[j]*x[ix+k];
          xvert[2*j+3] += bv[j]*x[ix+k+nbasis];
        }         

	//keep track of how many vertices are on right side of cut boundary
        printf("e %i, v %i %e %e \n",i,j,xvert[2*j],xvert[2*j+1]); 
        if(xvert[2*j]<x0) tally[j]=1;
      }

      // check if element has vertices on both sides of x=x0 then 
      // track which vertices will be cut and which kept
      sum=0;
      for(j=0;j<nfp;j++) sum+=tally[j];
      if(sum>0 && sum<nfp){
        l=0;
        m=0;
        for(j=0;j<nfp;j++){
          if(tally[j]==1){
            vcut[l++] = j;
          } else{
            vorig[n++] = j;
          }
        }
      }

      // need to store xcut (cut cell corner locations)
      // Currently only works for p=1 since that's all we need
      //
      // tri cut region (1 original node, 2 intersect pts)
      if(sum==1){ 
        xtmp[0]= xvert[2*vcut[0]];
        xtmp[0+3] = xvert[2*vcut[0]+3];

        // form edges between vcut and vorig vertices and get intersect pts
        // with cut boundary
        //
        // Find intercept by doing ycut = a*x0+b
 	
        // e1 = vcut - vorig[0]      
        a = (xvert[2*vorig[0]+3] - xvert[2*vcut[0]+3])/(xvert[2*vorig[0]] - xvert[2*vcut[0]]);
        b = xvert[2*vorig[0]+3]-a*xvert[2*vorig[0]]; 
        ycut1=a*x0+b; 
 
        xtmp[1]   = x0;
        xtmp[1+3] = ycut1; 
        
	// e2 = vcut - vorig[1]      
        a = (xvert[2*vorig[1]+3] - xvert[2*vcut[0]+3])/(xvert[2*vorig[1]] - xvert[2*vcut[0]]);
        b = xvert[2*vorig[1]+3]-a*xvert[2*vorig[1]]; 
        ycut2=a*x0+b; 

        xtmp[2]   = x0;
        xtmp[2+3] = ycut2; 

        // double check for positive jacobian XXX
        checkJac(xtmp,det); 
        if(det>0){
          xcut[m*3*2]   = xtmp[0];
          xcut[m*3*2+1] = xtmp[1];
          xcut[m*3*2+2] = xtmp[2];
          xcut[m*3*2+3] = xtmp[3];
          xcut[m*3*2+4] = xtmp[4];
          xcut[m*3*2+5] = xtmp[5];
        }
        else{
          xcut[m*3*2]   = xtmp[0];
          xcut[m*3*2+1] = xtmp[2];
          xcut[m*3*2+2] = xtmp[1];
          xcut[m*3*2+3] = xtmp[3];
          xcut[m*3*2+4] = xtmp[5];
          xcut[m*3*2+5] = xtmp[4];
        }
 
        cut2e[n]=i; //store id of orig elem
        n++; 
        necut++;
      }
      // quad cut region (2 original nodes, 2 intersect pts)
      // need to form 2 triangles here
      else if(sum==2){ 
	//First triangle (top cut node, 2 intersect nodes)
	if(xvert[2*vcut[0]+3]>xvert[2*vcut[1]+3]){
          xtmp[0]= xvert[2*vcut[0]];
          xtmp[0+3] = xvert[2*vcut[0]+3];
        }
        else{
          xtmp[0]= xvert[2*vcut[1]];
          xtmp[0+3] = xvert[2*vcut[1]+3];
        }

        // form edges between vcut and vorig vertices and get intersect pts
        // with cut boundary
        //
        // Find intercept by doing ycut = a*x0+b
 	
        // e1 = vcut - vorig[0]      
        a = (xvert[2*vorig[0]+3] - xvert[2*vcut[0]+3])/(xvert[2*vorig[0]] - xvert[2*vcut[0]]);
        b = xvert[2*vorig[0]+3]-a*xvert[2*vorig[0]]; 
        ycut1=a*x0+b; 
 
        xtmp[1]   = x0;
        xtmp[1+4] = ycut1; 
        
	// e2 = vcut - vorig[1]      
        a = (xvert[2*vorig[1]+3] - xvert[2*vcut[0]+3])/(xvert[2*vorig[1]] - xvert[2*vcut[0]]);
        b = xvert[2*vorig[1]+3]-a*xvert[2*vorig[1]]; 
        ycut2=a*x0+b; 

        xtmp[2]   = x0;
        xtmp[2+4] = ycut2; 

        // double check for positive jacobian XXX
        checkJac(xtmp,det);
        if(det>0){
          xcut[m*3*2]   = xtmp[0];
          xcut[m*3*2+1] = xtmp[1];
          xcut[m*3*2+2] = xtmp[2];
          xcut[m*3*2+3] = xtmp[3];
          xcut[m*3*2+4] = xtmp[4];
          xcut[m*3*2+5] = xtmp[5];
        }
        else{
          xcut[m*3*2]   = xtmp[0];
          xcut[m*3*2+1] = xtmp[2];
          xcut[m*3*2+2] = xtmp[1];
          xcut[m*3*2+3] = xtmp[3];
          xcut[m*3*2+4] = xtmp[5];
          xcut[m*3*2+5] = xtmp[4];
        }
 
        cut2e[n]=i; //store id of orig elem
        n++; 
        necut++;

	//Second triangle (both cut nodes, bottom intersect nodes)
        xtmp[0]= xvert[2*vcut[0]];
        xtmp[0+3] = xvert[2*vcut[0]+3];
        xtmp[1]= xvert[2*vcut[1]];
        xtmp[1+3] = xvert[2*vcut[1]+3];
        xtmp[2] = x0; 
	if(ycut1<ycut2){
          xtmp[2+3] = ycut1; 
        }
        else{
          xtmp[2+3] = ycut2; 
        }
	
        // double check for positive jacobian XXX
        checkJac(xtmp,det);
        if(det>0){
          xcut[m*3*2]   = xtmp[0];
          xcut[m*3*2+1] = xtmp[1];
          xcut[m*3*2+2] = xtmp[2];
          xcut[m*3*2+3] = xtmp[3];
          xcut[m*3*2+4] = xtmp[4];
          xcut[m*3*2+5] = xtmp[5];
        }
        else{
          xcut[m*3*2]   = xtmp[0];
          xcut[m*3*2+1] = xtmp[2];
          xcut[m*3*2+2] = xtmp[1];
          xcut[m*3*2+3] = xtmp[3];
          xcut[m*3*2+4] = xtmp[5];
          xcut[m*3*2+5] = xtmp[4];
        }
 
        cut2e[n]=i; //store id of orig elem
        n++; 
        necut++;

      }
  }


}
