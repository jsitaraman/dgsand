

void cross(double *a,double *b,double *c,int d)
{
   a[0]=b[1]*c[2]-b[2]*c[1];
   a[1]=b[2]*c[0]-b[0]*c[2];
   if (d == 3) { 
   a[2]=b[0]*c[1]-b[1]*c[0];
   }
   return;
}
  
void axb(double* a, double* x, double* b, int d)
{
  int i,j,ind; 
  for(i=0;i<d;i++){
    b[i]=0; 
    for(j=0;j<d;j++){
      ind = i*d+j; 
      b[i] += a[ind]*x[j];
    }
  }
}

double total_area(double *detJ, int etype, int p, int d, int nelem)
{
  int i,j,m;
  double totalArea=0;
  double wgt;  
  m=0;
  for(i=0;i<nelem;i++) 
    for(j=0;j<ngElem[etype][p];j++)
      {
	wgt=0.5*gauss[etype][p2g[etype][p]][(d+1)*j+2];
	totalArea+=(wgt*detJ[m]);
	m++;
      }
}

void CellCoordInterp(double *xsource, double *ijk, 
                     double *xdest, double *Jinvdest, double *rst,
                     int d, int e, int p, int ismerged)
//void CellCoordInterp(double *xdest, int d, int e, int p, double* Jinvdest, 
//                     double *ijk, double *xsource, 
//                     double* rst, int ismerged) 
//
// Routine that converts local coordinates from one cell to another
//
// Inputs: 
//   xdest = global coords of destination tri
//   d = number of spat dims
//   e = original element type
//   Jinvdest = destination element inv Jac
//   ijk = rst value of desired pt in source tri
//   xsource = global coords of source tri
//   ismerged = whether or not cell is merged to another
//
// Outputs
//   rst = rst value of desired point on destiation tri
{
 
  // notation: 
  //   xyz = global coord sys
  //   ijk = local coord sys of source tri
  //   rst = local coord sys of destination tri

  int b,w,i,j,l,ld,ij,m,n;
  int nbasisx=order2basis[e][1]; 
  int nbasis=order2basis[e][p]; 
  double jsource[d*d],mat[d][d],ijac[d][d],det;
  double xysource[d];
  double bd[nbasis][d]; 
  l=m=ld=ij=0;
  det = -1000; 

  //get xyz coord of nodes on the two elements
  //corresponding to rst = 0,0,1
  double x0[2] = {0.0, 0.0}; 
  double x0source[2] = {0.0, 0.0}; 
  double u[2] = {0,0};
  for(i=0;i<nbasisx;i++){
    x0source[0] = x0source[0] + xsource[i]*basis[e][i](u);
    x0source[1] = x0source[1] + xsource[i+nbasisx]*basis[e][i](u);
  }
  for(i=0;i<nbasis;i++){
    x0[0] = x0[0] + xdest[i]*basis[e][i](u);
    x0[1] = x0[1] + xdest[i+nbasis]*basis[e][i](u);
  }

  // Compute p=1 jacobian for subcell element, dx/di
  // only p=1 required for this regardless of element p 
  JacP1Tri(jsource,xsource,&det); 

  // Compute the global x y location of this source cell quad pt
  // [xysource] = [jac_source]*[ij] + x0_source
  axb(jsource,ijk,xysource,d);
  for(j=0;j<d;j++) xysource[j] = xysource[j] + x0source[j] - x0[j];
//  printf("\t\t\txy = %f %f\n",xysource[0],xysource[1]);
//  printf("\t\t\tJinvdest = %f %f; %f %f\n",Jinvdest[0],Jinvdest[1],Jinvdest[2],Jinvdest[3]); 

  // Compute rst coord of source cell quad
  // [rs] = [Jinv_dest][xysource - x0_dest]
  axb(Jinvdest,xysource,rst,d);
//  printf("\t\t\trs = %f %f\n",rst[0],rst[1]);

  if(ismerged==0 && (rst[0]<-1e-13 || rst[1]<-1e-13 || rst[0] > 1+1e-10 || rst[1] > 1+1e-10)){
    axb(jsource,ijk,xysource,d);
    printf("ERROR: \n\txyz loc = %.16e %.16e\n",xysource[0]+x0source[0],xysource[1]+x0source[1]); 
    printf("\tijk = %.16e %.16e\n\trst = %.16e %.16e\n",ijk[0],ijk[1],rst[0],rst[1]);

    printf("\tx0 = (%f %f),",x0[0],x0[1]); 
    u[0] = 1; 
    u[1] = 0;
    x0[0]=x0[1] = 0.0;  
    for(i=0;i<nbasis;i++){
      x0[0] = x0[0] + xdest[i]*basis[e][i](u);
      x0[1] = x0[1] + xdest[i+nbasis]*basis[e][i](u);
    }
    printf("\t(%f %f)",x0[0],x0[1]); 
    u[0] = 0; 
    u[1] = 1;
    x0[0]=x0[1] = 0.0;  
    for(i=0;i<nbasis;i++){
      x0[0] = x0[0] + xdest[i]*basis[e][i](u);
      x0[1] = x0[1] + xdest[i+nbasis]*basis[e][i](u);
    }
    printf("\t(%f %f)\n",x0[0],x0[1]); 

    printf("\tx0source = (%f %f),",x0source[0],x0source[1]); 
    u[0] = 1; 
    u[1] = 0;
    x0source[0]=x0source[1] = 0.0;  
    for(i=0;i<nbasisx;i++){
      x0source[0] = x0source[0] + xsource[i]*basis[e][i](u);
      x0source[1] = x0source[1] + xsource[i+nbasisx]*basis[e][i](u);
    }
    printf("\t(%f %f),",x0source[0],x0source[1]); 
    u[0] = 0; 
    u[1] = 1;
    x0source[0]=x0source[1] = 0.0;  
    for(i=0;i<nbasisx;i++){
      x0source[0] = x0source[0] + xsource[i]*basis[e][i](u);
      x0source[1] = x0source[1] + xsource[i+nbasisx]*basis[e][i](u);
    }
    printf("\t(%f %f)\n",x0source[0],x0source[1]); 

    printf("\torig Jinvdest = %f %f ; %f %f\n",Jinvdest[0],Jinvdest[1],Jinvdest[2],Jinvdest[3]);     
    exit(1);  
  }
}

void BasesVCut(double *x, double *Jinv,double *detJ,
               double *xcut, double *bvcut, double *bvdcut, 
               double *Jinvcut, double *detJcut,
               int d, int e, int p,int icut,int eid,int pid,int ismerged)

// Routine that builds up the bases and derivs of the interior cut cell quad pts
// Note:
// 	x = global coord of PARENT element
// 	xcut = global coord of cut cell
{
  int g = p2g[e][p];
  int l,n,m,ld,ii,w,b,i,j; 
  int nbasis1=order2basis[e][1];       
  int nbasis=order2basis[e][p+(p==0)];
  double u[2],bd[nbasis][d],mat[d][d],jac[d][d],det;

  //Loop over vol quad pts  
  double ijk[d]; 
  l=n=m=ld=ii=0;

  // debug
  if(ismerged) printf("\ndebug BasesVCut: eid, pid, icut, ismerged = %i %i %i %i\n",eid,pid,icut,ismerged);
  for(w=0; w<ngElem[e][p]; w++){
    // get local coord of quad pt on cut cell 
    for(j=0; j<d;j++) ijk[j] = gauss[e][g][(d+1)*w+j]; //double check inputs to gauss XXX
    
    // get bases and derivs at quad pt 
    CellCoordInterp(xcut, ijk, x, Jinv, u, d, e, p,ismerged); 
    //CellCoordInterp(x,d,e,p,Jinv,ijk,xcut,u,ismerged); 

    //store bvcut at each quad pt
    for(b=0;b<nbasis;b++){ // loop over bases
      for(j=0;j<d;j++){
	    bd[b][j]=basis_d[e][b*d+j](u); // accumulate bases derivs at quad pt
      }
      if (p > 0) bvcut[l++]=basis[e][b](u); // filled in as bv[nGL][nbasis]
    }
    if (p==0) bvcut[l++]=1;

    if(ismerged){
      for(b=0;b<nbasis;b++)
        for(j=0;j<d;j++)
	  printf("\tbd[%i][%i] = %f\n",b,j,bd[b][j]);
    }

    //build jacobian dx/dr of cut cell
    //How to get this for p>1? Only have 3 xcut cells
    for(i=0;i<d;i++){
      for(j=0;j<d;j++){
	mat[i][j]=xcut[i*nbasis1]*bd[0][j];
	for(b=1;b<nbasis1;b++)
  	  mat[i][j]+=xcut[i*nbasis1+b]*bd[b][j];
      }
    }

    //invert jacobian (stored in jac) and get detJ
    if (d==2) invmat2x2(mat,jac,det);
    for(i=0;i<d;i++)
      for(j=0;j<d;j++){ 
        Jinvcut[ii++]=jac[i][j];
      }

    // get basis derivs dN/dx = dN/dr*(dx/dr)^-1
    // Use full cell inv Jacobiain b/c we're using full cell shape functions
    if (p > 0) {
      for(b=0;b<nbasis;b++)
	for(i=0;i<d;i++){
	   bvdcut[ld]=Jinv[i]*bd[b][0];
	   for(j=1;j<d;j++)
	     bvdcut[ld]+=Jinv[j*d+i]*bd[b][j]; 
	   ld++;
        }
    }
    else {
      for(i=0;i<d;i++) bvdcut[ld++]=0;
    }
    detJcut[n++]=det;
  }
}

void CutFaceWeights(double *x, double *Jinv, int pc, int* iptr, double *xcut, double *bfcutL, double *bfdcutL, double *bfcutR, double *bfdcutR, double *Jinvcut, double *faceWeight, int d, int e, int p, int* cut2neigh, int* cutoverset, int iorig, int icut, int* elemParent)
{
  int b,w,i,j,ij,l,ld,n,m,f,f1,pid,neigh, ismerged;
  double uL[d],uR[d],v,wgt;
  double mat[d][d],jacL[d][d],jacR[d][d],det;
  int nbasis=order2basis[e][p+(p==0)];
  double bdL[nbasis][d];  // basis derivative  
  double bdR[nbasis][d];  // basis derivative  
  int nfaces=facePerElem[e];
  int g=p2gf[e][p]; // gauss quadrature type for this element type
  double xx,yy;
  // ok, we are not going above 3D now
  double Ja[3];
  double Jb[3]={0,0,1};

  // for every face
  double ijk[2];


  printf("debug: cutfaceweights: icut = %i, iorig = %i, iparent = %i\n",icut,iorig,elemParent[iorig]);
  m=l=0;
  ld=ij=0;
  for(f=0;f<nfaces;f++){    
    // for every Gauss-point on this sub face
    for(w=0;w<ngGL[e][p];w++){
      printf("\t face %i weight %i\n",f,w);
      v=gaussgl[e][g][(d)*w];  // gauss location from 0 to 1
	    wgt=gaussgl[e][g][(d)*w+1];	 // gauss weight

      // this is specific to 2-D 
	    // remember fill order will be reversed for element
      // sharing this edge   
      f1=(f+1)%nfaces;

    	// do 1D interp to get rst coord (on sub element) of quad pt
	    ijk[0]=(1-v)*eloc[e][p][d*f]  +v*eloc[e][p][d*f1]; 
    	ijk[1]=(1-v)*eloc[e][p][d*f+1]+v*eloc[e][p][d*f1+1];

      // convert ijk coord on sub elem to rst coord on full L side elem
      pid = elemParent[iorig]; 
      ismerged = 0; 
      if(iorig!=elemParent[iorig]) ismerged = 1; 
      CellCoordInterp(xcut, ijk, x+iptr[pc*pid+1], Jinv+iptr[pc*pid+8],
                      uL, d, e, p, ismerged); 
//      CellCoordInterp(x+iptr[pc*pid+1],d,e,p,Jinv+iptr[pc*pid+8],ijk,xcut,uL,ismerged); 
      printf("\t\tuL = %f %f\n",uL[0],uL[1]);

      // convert ijk coord on sub elem to rst coord on full R side elem
      neigh = cut2neigh[f];
      ismerged = 0; 
      if(neigh!=elemParent[neigh]) ismerged = 1; 
      printf("\t\tneigh = %i, neighParent = %i\n",neigh,elemParent[neigh]);
      if(neigh>-1){
    	  pid = elemParent[neigh];
        CellCoordInterp(xcut, ijk, x+iptr[pc*pid+1], Jinv+iptr[pc*pid+8],
                        uR, d, e, p, ismerged); 
        //CellCoordInterp(x+iptr[pc*pid+1],d,e,p,Jinv+iptr[pc*pid+8],ijk,xcut,uR,ismerged); 
        printf("\t\tuR = %f %f\n",uR[0],uR[1]);
	     }

  	   // get the bases and derivs at gauss legendre pts
	     for(b=0;b<nbasis;b++){
    	   for(j=0;j<d;j++){
        	 bdL[b][j]=basis_d[e][b*d+j](uL); 
           if(neigh>=-1)
             bdR[b][j]=basis_d[e][b*d+j](uR); 
         }

  	     // basis[][](u) computes basis function at u
    	   // Skip this step if the face is an overset boundary XXX
	       if (p > 0){
          	bfcutL[l]=basis[e][b](uL);  // filled as bf[nfaces][nGL][nbasis]
          	if(neigh>-1){
      	      bfcutR[l]=basis[e][b](uR);  // filled as bf[nfaces][nGL][nbasis]
            }
          	l++;
	       } 
	     }
    	 if (p==0){
	       bfcutL[l]=1.0;
    	   bfcutR[l]=1.0;
	       l++; 
    	 }

    	 // build the jacobians of the cut cell
	     for(i=0;i<d;i++){	    
         faceWeight[d*m+i]=0;
    	   Ja[i]=0;
	       for(j=0;j<d;j++){
	         mat[i][j]=xcut[i*3]*bdL[0][j]; 
  	       for(b=1;b<3;b++)
	           mat[i][j]+=xcut[i*3+b]*bdL[b][j]; 
	       } 

    	   // need a cross product here 
         for(j=0;j<d;j++){	       
           Ja[i]+=(mat[i][j]*face2elem[e][d*f+j]); // face2elem is...?
		       //TODO add Jb[i] calculation for 3D elements here
     	   }
       }

    	 // do faceWeight = Ja x zhat
	     // This gives me the normal vector
    	 cross(&(faceWeight[2*m]),Ja,Jb,d); 

       printf("\tmat = %f %f %f %f\n",mat[0][0],mat[0][1],mat[1][0],mat[1][1]);
       printf("\tJa = %f %f %f\n",Ja[0],Ja[1],Ja[2]); 
       printf("\tfwcut[%i-%i] = %f %f \n",2*m,2*m+1,faceWeight[2*m],faceWeight[2*m+1]); 

       // get [dx/dr]^-1
       if (d==2) invmat2x2(mat,jacL,det);

//Commenting all this out b/c currently not using viscous
/*
 	    // build the jacobians of the R neighbor cell 
//   	    How do I account for the wrong arc length b/c we're using x and not xcut

  	    if(eid!=-1){
  	      for(i=0;i<d;i++){	    
	        for(j=0;j<d;j++){
                  mat[i][j] =  x[iptr[eid*pc+1]+i*3]*bdR[0][j]; 
	          for(b=1;b<3;b++)
	            mat[i][j]+=x[iptr[eid*pc+1]+i*3+b]*bdR[b][j];  
  	        } 
              }
              if (d==2) invmat2x2(mat,jacR,det);

              // XXX need to multiply det by cut length? XXX

              // get [dx/dr]^-1
              if (d==2) invmat2x2(mat,jacR,det);
	    }

            // compute [dr/dx][dN/dr]
	    if (p > 0) {
 	      for(b=0;b<nbasis;b++){ 
		for(i=0;i<d;i++){

		  bfdcutL[ld]=jacL[0][i]*bdL[b][0];
		  for(j=1;j<d;j++)
		    bfdcutL[ld]+=jacL[j][i]*bdL[b][j];

		  if(eid!=-1){
  		    bfdcutR[ld]=jacL[0][i]*bdR[b][0];
		    for(j=1;j<d;j++)
		      bfdcutR[ld]+=jacL[j][i]*bdR[b][j]; // what to do abou jac here? should it be jacL or R?
		  }
		  else{
  		    bfdcutR[ld] = 0.0; 
		  }

		  ld++;
		}
	      }
	    }
	    else {
              for(i=0;i<d;i++){
		 bfdcutL[ld]=0;
		 bfdcutR[ld]=0; 
		 ld++;
	      }
	    }
*/
// MEMORY LEAK HERE, fix later
//	    for(i=0;i<d;i++)
//              for(j=0;j<d;j++)	     
//  	        Jinvcut[ij++]=jacL[i][j];
	    m++;
	}// loop over gauss pts
//     } // if cut overset face
   } //loop over faces
}

void cut_mass_matrix(double *M, double *x, double *Jinv, double *xcut, double *detJcut, int d, int e, int p, int ismerged)
{
  int i,j,ij,w,b,ii,jj;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];
  double mat[d][d],jac[d][d];
  double ijk[d],u[d],wgt,det;
  int g=p2g[e][p+1]; // gauss point type for this element type, use one
                     // order higher making sure mass matrix is exact

  if (p > 0 ){

    // Get jacobian info for the original element at p+1 quad pts
    // This is needed to interp between cut cell coords and parent elem coords
    for(i=0;i<nbasis;i++)
      for(j=0;j<nbasis;j++)
	{
	  ij=nbasis*i+j;
	  for(w=0;w<ngElem[e][p+1];w++) // 
	    {
	      for(jj=0;jj<d;jj++)
      		ijk[jj]=gauss[e][g][(d+1)*w+jj];

	      wgt=gauss[e][g][(d+1)*w+2]*detJcut[0];

	      // could use bv here instead of reevaluating 
	      // revaluate this if mesh is deforming 

        // subtract cut region from original mass matrix
        CellCoordInterp(xcut, ijk, x, Jinv, u, d, e, p, ismerged); 
        //CellCoordInterp(x,d,e,p,Jinv,ijk,xcut,u,ismerged); //,invJcut,detJcut); 
	      M[ij]-=(wgt*basis[e][i](u)*basis[e][j](u));
	    }
  	}
  }
  /// XXX What to do here? 
  else {
    // p=0 
      M[0]-=0.5*detJcut[0];
  }
}

void mass_matrix(double *M, double *x, double* Jinv, double* detJ, int d, int e, int p, int debug, int eid, int pid, int* iptr, int pc)
{
  int i,j,ij,w,b,ii,jj;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];
  double mat[d][d],jac[d][d];
  double u[d],utmp[d],wgt,det;
  int g=p2g[e][p+1]; // gauss point type for this element type, use one
                     // order higher making sure mass matrix is exact

  double* xparent = x+iptr[pc*pid+1]; 
  double* Jinvparent = Jinv+iptr[pc*pid+4];
  double* xcur = x+iptr[pc*eid+1]; 
  int ismerged = (pid!=eid);
  if (p > 0 ) {
    for(i=0;i<nbasis;i++)
      for(j=0;j<nbasis;j++){
        for(w=0;w<ngElem[e][p+1];w++){
          // local coord of quad pt in current element
          for(jj=0;jj<d;jj++)
            utmp[jj]=gauss[e][g][(d+1)*w+jj]; 

            // convert to local coord of parent element
            // (parent != self for merged elems)
            if(ismerged){
              CellCoordInterp(xcur, utmp, xparent, Jinvparent, 
                              u, d, e, p, ismerged); 
              //CellCoordInterp(xparent,d,e,p,Jinvparent,utmp,xcur,u,ismerged);
              if(debug){
                double x0cur[2]={0,0},x0par[2]={0,0};
                int nbasisx=order2basis[e][1]; 
                double v[2] = {0,0};
                for(int k=0;k<nbasisx;k++){
                  x0par[0] += xparent[k]*basis[e][k](v);
                  x0par[1] += xparent[k+nbasisx]*basis[e][k](v);
                  x0cur[0] += xcur[k]*basis[e][k](v);
                  x0cur[1] += xcur[k+nbasisx]*basis[e][k](v);
                }
                printf("\t\t\tJinvparr = %f %f %f %f\n",Jinvparent[0],Jinvparent[1],Jinvparent[2],Jinvparent[3]);
                printf("\t\t\tx0par = %f %f\n\t\t\tx0cur = %f %f\n",
                        x0par[0],x0par[1],x0cur[0],x0cur[1]);            
              }
            } 
	          else{
              for(jj=0;jj<d;jj++) u[jj]=utmp[jj];
            }

            // get weight of current element
            // and update parent mass mastrix
            wgt=gauss[e][g][(d+1)*w+2]*detJ[iptr[pc*eid+5]];
            M[iptr[pc*pid+10]+nbasis*i+j]+=(wgt*basis[e][i](u)*basis[e][j](u));
        } // loop quad pts
      } // loop j shp functions
  }
  else {
    /* p=0 */
      for(b=0;b<nbasis;b++)
       for(jj=0;jj<d;jj++)
          bd[b][jj]=basis_d[e][b*d+jj](u);
      
      for(ii=0;ii<d;ii++)
	{
	  for(jj=0;jj<d;jj++)
	    {
	      mat[ii][jj]=x[ii*nbasis]*bd[0][jj];
	      for(b=1;b<nbasis;b++)
		mat[ii][jj]+=x[ii*nbasis+b]*bd[b][jj];
	    }
	}
      if (d==2) invmat2x2(mat,jac,det);
      M[0]=0.5*det;
  }
}

void Jacobian(double *x,double *Jinv, double * detJ,
              int d, int e, int p)
{
  int b,w,i,j,l,ld,ij,m,n;
  double u[d];
  double mat[d][d],jac[d][d],det;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];  // basis derivative  
  int g=p2g[e][p];  // gauss quadrature data for this element type

  l=n=m=ld=ij=0;
  for(w=0;w<ngElem[e][p];w++) // loop over gauss pts
    {

      for(j=0;j<d;j++) // get rs coordinates of gauss pt
        u[j]=gauss[e][g][(d+1)*w+j];

      // accumulate basis derivs at quad pt
      // in local element system
      for(b=0;b<nbasis;b++)
        for(j=0;j<d;j++)
          bd[b][j]=basis_d[e][b*d+j](u);
     
      //build jacobian dx/dr for local element
      for(i=0;i<d;i++){
	      for(j=0;j<d;j++){
          mat[i][j]=x[i*nbasis]*bd[0][j];
          for(b=1;b<nbasis+(nbasis==1);b++)
          	mat[i][j]+=x[i*nbasis+b]*bd[b][j];
        }
      }

      for(i=0;i<d;i++)
        for(j=0;j<d;j++)
          //invert jacobian (stored in jac) and get detJ
          if (d==2) invmat2x2(mat,jac,det);
          for(i=0;i<d;i++)
            for(j=0;j<d;j++) 
              Jinv[ij++]=jac[i][j];
      detJ[n++]=det;
    }
} 

void VolWeights(double *x,double *bv, double *bvd, double *Jinvparent, 
              double *detJ, int* iptr, int pc, int eid, int pid, 
              int d, int e, int p)
{
  int b,w,i,j,l,ld,ij,m,n;
  double ucur[d],upar[d];
  double mat[d][d],jac[d][d],det;
  int nbasis=order2basis[e][p+(p==0)];
  int g=p2g[e][p];  // gauss quadrature data for this element type

  double* xcur = x+iptr[pc*eid+1]; 
  double* xparent = x+iptr[pc*pid+1]; 
  double JinvT[d*d], bd[d];

  l=n=m=ld=ij=0;
  for(w=0;w<ngElem[e][p];w++){
    // get rs coordinates of gauss pt in local element
    for(j=0;j<d;j++) 
      ucur[j]=gauss[e][g][(d+1)*w+j];

    // convert local rst coordinates 
    // into parent elem system 
    if(pid!=eid){ 
      CellCoordInterp(xcur, ucur, xparent, Jinvparent, 
                      upar, d, e, p, (eid!=pid));
      //CellCoordInterp(xparent,d,e,p,Jinvparent,ucur,xcur,upar,(pid!=eid));
    }
    else{
      for(j=0;j<d;j++) upar[j] = ucur[j];
    }

    if (p > 0) {
      for(b=0;b<nbasis;b++){
        // fill in shape function value bv[nGL][nbasis]
        bv[l++]=basis[e][b](upar);

        // get bases and derivs in parent system 
        // dN/dx = trans(Jinv)*[dN/dr;dN/ds]
        for(i=0;i<d;i++){
          bd[i] = basis_d[e][b*d+i](upar);           
          for(j=0;j<d;j++) JinvT[i*d+j] = Jinvparent[j*d+i];
        }
        axb(JinvT,bd,bvd+ld,d); 
        ld += 2; // increment counter up to next shp function b
      }
    }
    else {
      bv[l++]=1;
      for(i=0;i<d;i++) bvd[ld++]=0;
    }
  }

/*
  l=n=m=ld=ij=0;
  for(w=0;w<ngElem[e][p];w++){
    for(b=0;b<nbasis;b++){
      printf("w=%i,b=%i,bv = %f\n",w,b,bv[l++]);
      for(i=0;i<2;i++){
        printf("dN/dx%i = %f\n", i,bvd[ld++]); 
      }
    }
  }
*/

}

void FaceWeights(double *x, double *bf, double *bfd, double *JinvV, 
                 double *JinvF, double *faceWeight, 
              	 int d, int e, int p, int eid, int pid,
                 int* iptr, int pc)
{
  int b,w,i,j,ij,l,ld,n,m,f,f1;
  double ucur[d],upar[d],v,wgt;
  double mat[d][d],jac[d][d],det;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[d];  // basis derivative  
  int nfaces=facePerElem[e];
  int g=p2gf[e][p]; // gauss quadrature type for this element type
  double xx,yy;
  // ok, we are not going above 3D now
  double JinvT[d*d];
  double Ja[3];
  double Jb[3]={0,0,1};

  // for every face
  m=l=0;
  ld=ij=0;

  double* xcur    = x+iptr[pc*eid+1]; 
  double* xparent = x+iptr[pc*pid+1]; 
  double* Jinvcur = JinvV+iptr[pc*eid+4];
  double* Jinvpar = JinvV+iptr[pc*pid+4];
  //printf("---------------\n");
  for(f=0;f<nfaces;f++)
    {
      // for every Gauss-point on this face
      for(w=0;w<ngGL[e][p];w++)
        {

          v=gaussgl[e][g][(d)*w];  // gauss location
          wgt=gaussgl[e][g][(d)*w+1];	 // gauss weight

          /* this is specific to 2-D */
          /* remember fill order will be reversed for element
             sharing this edge */  
          f1=(f+1)%nfaces;

          // get location of gauss pt in local and parent element coord sys
          ucur[0]=(1-v)*eloc[e][p][d*f]  +v*eloc[e][p][d*f1];
          ucur[1]=(1-v)*eloc[e][p][d*f+1]+v*eloc[e][p][d*f1+1];
          if(eid!=pid){ // cell merged
            // convert ucur into parent element coord rst
            CellCoordInterp(xcur, ucur, xparent, Jinvpar, 
                            upar, d, e, p, (eid!=pid));
            //CellCoordInterp(xparent,d,e,p,Jinvpar,ucur,xcur,upar,(pid!=eid));
          }
          else{ // not cell merged
            for(i=0;i<d;i++) upar[i] = ucur[i];
          }
	  
printf("\tf=%i, w=%i, ucur = %f %f, upar = %f %f\n",f,w,ucur[0],ucur[1],upar[0],upar[1]); 

          // Assuming straight edge elements
          // Jacobian is constant everywhere in element
          // Using JinvV that's already precalculated
          for(i=0;i<d;i++)
            for(j=0;j<d;j++){
              JinvF[ij++]=Jinvpar[d*i+j];
              mat[i][j] =Jinvcur[d*i+j];
            }
          invmat2x2(mat,jac,det); // get Jac of current elem

          // Get shape function value at gauss pt 
          // using parent bases
          for(b=0;b<nbasis;b++){
            if (p > 0) bf[l++]=basis[e][b](upar);  // filled as bf[nfaces][nGL][nbasis]
printf("\tf=%i, w=%i, b=%i, bf[%i] = %f\n",f,w,b,l-1,bf[l-1]);
          }

          // Get shape function derivs at gauss pt
          // using parent bases
          // dN/dx = trans(Jinv)*[dN/dr;dN/ds]
          for(b=0;b<nbasis;b++){
            for(i=0;i<d;i++){
              bd[i] = basis_d[e][b*d+i](upar);           
              for(j=0;j<d;j++) JinvT[i*d+j] = Jinvpar[j*d+i];
            }
            axb(JinvT,bd,bfd+ld,d); 
printf("\tf=%i, w=%i, b=%i, dir=%i, bfd[%i, %i] = %f %f\n",f,w,b,i,ld,ld+1,bfd[ld],bfd[ld+1]);
            ld += 2; // increment counter up to next shp function b
          }
  
          // Get edge vector in physical space
          for(i=0;i<d;i++){
            Ja[i]=0;
            for(j=0;j<d;j++)
              // face2elem is edge vector in rst of current element
              Ja[i]+=(jac[i][j]*face2elem[e][d*f+j]); 
          }
printf("\tf=%i, w=%i, Jinv = %f %f %f %f\n",f,w,Jinvcur[0],Jinvcur[1],Jinvcur[2],Jinvcur[3]);
          // do faceWeight = Ja x zhat
          // This gives me the normal vector in physical space
          cross(&(faceWeight[2*m]),Ja,Jb,d); 
printf("\tf=%i, w=%i, Ja = %f %f\n",f,w,Ja[0],Ja[1]);
printf("\tf=%i, w=%i, faceWeight = %f %f\n",f,w,faceWeight[2*m],faceWeight[2*m+1]);
          
          // 
          m++;
        } // loop over gauss
    } // loop over faces
}

void COMPUTE_GRID_METRICS(double *x, double *bv, double *bvd,double *JinvV, 
			  double *detJ,double *bf, double *bfd, double *JinvF, double *faceWeight,
			  int *iptr, int *elemParent, int d, int e, int p, int nelem, 
        int pc, int* iblank, int imesh)
{
  int i,j,b,pid;
  int ip,ix,ibv,ibvd,ibf,ij,idetj,ibfd,ijf,ifw;
  
  // Need to compute jacobians for all full elements first
  for(i=0;i<nelem;i++)
    {
      iblank[i] = 0; 

      ip=pc*i;
      ix   =iptr[ip+1];
      ij   =iptr[ip+4];
      idetj=iptr[ip+5];

      printf("Debug Jacobian Elem %i, ij = %i\n",i,ij);
      Jacobian(x+ix,JinvV+ij,detJ+idetj,d,e,p); // basis on vol
    }
}

void COMPUTE_SHAPE(double *x, double *bv, double *bvd,double *JinvV, 
			  double *detJ,double *bf, double *bfd, double *JinvF, double *faceWeight,
			  int *iptr, int *elemParent, int d, int e, int p, int nelem, 
        int pc, int* iblank, int imesh)
{
  int i,j,b,pid;
  int ip,ix,ibv,ibvd,ibf,ij,idetj,ibfd,ijf,ifw,ijp;
  
  int nbasis=order2basis[e][p+(p==0)];
  // Now compute shape functions and derivs
  for(i=0;i<nelem;i++)
    {
      pid = elemParent[i];  // fully zero b/c haven't computed cut cells yet

printf("COMPUTESHAPE: ELEM %i, Parent %i\n",i,pid); 

      ip=pc*i;
      ibv  =iptr[ip+2];
      ibvd =iptr[ip+3];
      ibf  =iptr[ip+6];
      ibfd =iptr[ip+7];
      ijf  =iptr[ip+8];
      ifw  =iptr[ip+9];

      ijp =iptr[pc*pid+4]; // parent inverse jacobian
printf("\t ip %i, ibv %i, ibvd %i, ibf %i, ibfd %i, ijf %i, ifw %i\n",ip,ibv,ibvd, ibf, ibfd, ijf, ifw); 

      VolWeights(x, bv+ibv, bvd+ibvd,JinvV+ijp,detJ+idetj,iptr,pc,i,pid,d,e,p); // basis on vol

      FaceWeights(x,bf+ibf,bfd+ibfd,JinvV,JinvF+ijf,faceWeight+ifw,
                  d,e,p,i,pid,iptr,pc); // basis on face
    }
}

void COMPUTE_CUT_METRICS(double *x, double *JinvV, 
			 double *detJ,double *JinvF,
			 int *iptr, int d, int e, int p, int pc, int pccut,
			 double *xcut, double *bvcut, double *bvdcut, 
			 double *JinvVcut,double *JinvFcut,
			 double *detJcut, double *bfcutL, double *bfdcutL,  
			 double *bfcutR, double *bfdcutR,  
			 double *fwcut, int* iptrc, int necut, int* cut2e,
			 int* cut2neigh, int imesh, int* cutoverset, int* elemParent)
{
  int i,j,b,eid,pid;
  int ip,ix,ij,idetj,ijf;
  int cip,cix,cibv,cibvd,cibf,cij,cidetj,cibfd,cijf,cifw,cc2n,cico;

  // Metrics for the cut regions
  for(i=0;i<necut;i++)
  {
    eid = cut2e[i]; // Get elem ID of original elem 
    pid = elemParent[eid]; // parent element ID
    ip = pid*pc; 

    // Quantities for the parent element
    ix   =iptr[ip+1]; 
    ij   =iptr[ip+4];
    idetj=iptr[ip+5];

    ijf  =iptr[ip+8];

    // Quantities for the cut element
    cip = pccut*i;

    cix   =iptrc[cip+1]; // x coord
    cibv  =iptrc[cip+2]; // shp func      
    cibvd =iptrc[cip+3]; // shp func deriv
    cij   =iptrc[cip+4]; // Jinv
    cidetj=iptrc[cip+5]; // detJ

    cibf  =iptrc[cip+6]; // face shp function
    cibfd =iptrc[cip+7]; // face shp deriv 
    cijf  =iptrc[cip+8]; // JinvF
    cifw  =iptrc[cip+9]; // faceWeight
    cc2n  =iptrc[cip+12]; // cut2neigh  
    cico = iptrc[cip+12];

    // get bases at cut vol quad pts
    int ismerged=0; 
    if(eid!=pid) ismerged=1; 
    BasesVCut(x+ix, JinvV+ij, detJ+idetj, 
              xcut+cix, bvcut+cibv, bvdcut+cibvd,
              JinvVcut+cij, detJcut+cidetj,
              d,e,p,i,eid,pid,ismerged); 


    // basis on face
    printf("cut face weight cut elem %i\n",i);
    CutFaceWeights(x, JinvF,pc,iptr,xcut+cix,
  		   bfcutL+cibf,bfdcutL+cibfd,
  		   bfcutR+cibf,bfdcutR+cibfd,
                   JinvFcut+cijf,fwcut+cifw,d,e,p,
		   cut2neigh+cc2n,cutoverset+cico,eid,i,elemParent); 
  }

}

void MASS_MATRIX(double *mass,double *x, int *iptr, int d, int e, int p, int nelem, int pc, int imesh,
                 int* elemParent, double* Jinv, double* detJ)
{
  int i,pid,m;
  int ix,im,ixc;
  //compute mass matrix of all elements
  int nbasis=order2basis[e][p]; // XXX double check last argument, should be p=1 for cut cell
  int debug;

  // zero mass matrix
  for(i=0;i<nelem;i++){
    im=iptr[pc*i+10]; // parent mass matrix
    m = 0; 
    for(int j=0;j<nbasis;j++)
      for(int k=0;k<nbasis;k++){
        mass[im+m] = 0.0; 
        m++; 
      }
  }

  for(i=0;i<nelem;i++){
    pid=elemParent[i];
    debug = 1;
    if(debug) printf("\n\nDEBUG: Mesh %i, cell %i, parent %i:\n",imesh,i,pid); 
    ix=iptr[pc*i+1]; // current element x
    im=iptr[pc*pid+10]; // parent mass matrix
    mass_matrix(mass,x,Jinv,detJ,d,e,p,debug,i,pid,iptr,pc); 
  }

  /*// DEBUG 
  for(i=0;i<nelem;i++){
    pid=elemParent[i];
    if(pid!=i) debug = 1;
    if(debug){
        printf("\n");
        im=iptr[pc*i+10]; // parent mass matrix
        m = 0; 
        for(int j=0;j<nbasis;j++)
          for(int k=0;k<nbasis;k++){
            printf("Mesh %i, cell %i, full mass[%i] = %f\n",imesh,i,im+m,mass[im+m]);
            m++;
          }
    }
  }
  */
}

void CUT_MASS_MATRIX(double *mass,double *x, double *Jinv, int *iptr, double *xcut, double *detJcut, 
                     int *iptrc, int d, int e, int p, int nelem, int pc, int pccut, int necut, 
                     int* cut2e, int imesh, int* elemParent)
{
  int i,j,k,ld,eid,pid;
  int ix,im,ij,ixc,id;
  int ismerged;
  int nbasis=order2basis[e][p]; // XXX double check last argument, should be p=1 for cut cell
  //subtract out cut portion
  for(i=0;i<necut;i++){
      eid = cut2e[i]; // find original element id 
      pid = elemParent[eid]; // find original element id 
      ix=iptr[pc*pid+1];
      ij=iptr[pc*pid+4];
      im=iptr[pc*pid+10];

      ixc=iptrc[pccut*i+1];
      id=iptrc[pccut*i+5]; 

      ismerged=0;
      if(pid!=eid) ismerged=1; 
      cut_mass_matrix(&(mass[im]),&(x[ix]),&(Jinv[ij]),&(xcut[ixc]),&(detJcut[id]),d,e,p,ismerged);

      int debug = 1; 
      if(debug){
        int m = 0; 
        for(int j=0;j<nbasis;j++)
          for(int k=0;k<nbasis;k++){
            printf("Mesh %i, cell %i, cut mass[%i] = %f\n",imesh,i,im+m,mass[im+m]);
            m++;
          }
          printf("Mesh %i, cell %i, cut %i detJ = %f\n",imesh,eid,i,detJcut[id]); 
      }

  }
}

void findCentroid(double* x, double* centroid, int nbasis, int npf, int e, int d, int p)
{
  double bv,u[2];

  centroid[0] = 0.0; 
  centroid[1] = 0.0; 
  // get bases at rst = [0 0 ; 1 0; 0 1]
  for(int i=0;i<npf;i++){ // loop over vertices
    u[0]=eloc[e][1][d*i];
    u[1]=eloc[e][1][d*i+1];

    // get the vertex global coordinates of the original trying
    // by doing x = sum(Ni * xi)
    for(int j=0;j<nbasis;j++){ // accumulate bases and vertex coords
      bv = basis[e][j](u);
      centroid[0] += bv*x[j]/npf;
      centroid[1] += bv*x[j+nbasis]/npf;
    }         
  }
}

void FIND_PARENTS(double* x, int* iptr, int* elem2face, int* faces,
                  int* iblank, int* cellmerge, int* elemParent,
                  int d, int e, int p, int nelem, int pc, int imesh)
{
  double centroid[2],ncentroid[2],u[2]; 
  double dist,ndist,ndist2; 
  int npf = facePerElem[e];
  int nbasis=order2basis[e][p];
  int ix;
  int keep,neigh,cneigh,eid,fid;
  int maxiter,iter;

  printf("\n===================\n");
  printf("Mesh %i Parents\n",imesh);
  printf("\n===================\n");

  for(int i=0;i<nelem;i++){
    if(cellmerge[i]==2){ // severely cut element that needs merging
      elemParent[i] = -1; 
      cneigh = -1; 

      ix=iptr[pc*i+1];
      findCentroid(x+ix,centroid,nbasis,npf,e,d,p);

      // loop face neighbors
      maxiter = 2; // max number of levels of face neighbors to check
      iter = 1;     
      eid = i; // start with current element
      while(elemParent[i]==-1 && iter<=maxiter){ 
        keep = -1; 
        ndist = 1e16; 
        ndist2 = ndist;
        for(int j=0;j<npf;j++){
          //get faceID and neighbor elem id
          fid = abs(elem2face[3*eid+j])-1;
          neigh = faces[6*fid+2] == eid ? faces[6*fid+4] : faces[6*fid+2];

          if(neigh>=0 && neigh != i){
            //get neighbor elemID   
            findCentroid(x+iptr[pc*neigh+1],ncentroid,nbasis,npf,e,d,p);
            dist =  (ncentroid[0]-centroid[0])*(ncentroid[0]-centroid[0]);
            dist += (ncentroid[1]-centroid[1])*(ncentroid[1]-centroid[1]);

            //find nearest unblanked, stable element 
            if(cellmerge[neigh]!=2 && dist<ndist && iblank[neigh]!=1){
       	      ndist = dist; 
       	      keep = neigh;
            }

      	    // save closest unblanked neighbor in case this doesn't work
      	    if(iblank[neigh]!=1 && dist<ndist2){
      	      ndist2 = dist;
      	      cneigh = neigh;
            }
      	  } 
        }// loop face neighbors

      	// check to see if this level worked
      	if(keep != -1) elemParent[i] = keep;

        // if not, update eid and try again 
        eid = cneigh;
        iter++;
      } // loop through levels of face neighbors

      // if all face neighbors are STILL cut severely or blanked, give up
      if(elemParent[i]==-1){ 
      	printf("\nERROR: Bad hole cutting. Cannot find good parent element. Exiting.\n");
        exit(1); 
      }
    }
    else{ // elem is not severely cut
      elemParent[i] = i; 
    }

    //Debug outputs
    printf("Element %i, Parent is %i\n",i,elemParent[i]);
  } // loop nelem
}

