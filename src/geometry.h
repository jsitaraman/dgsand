

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
	wgt=gauss[etype][p2g[etype][p]][(d+1)*j+2];
	totalArea+=(wgt*detJ[m]);
	m++;
      }
}

double total_area_cut(double *detJ, int *iptrc, int pccut, int etype, int p, int d, int necut)
{
  int i,j,m,id;
  double totalArea=0;
  double wgt;
  m=0;
  for(i=0;i<necut;i++)
    {
      id=iptrc[pccut*i+5];
      for(j=0;j<ngElem[etype][p];j++)
	{
	  wgt=gauss[etype][p2g[etype][p]][(d+1)*j+2];
	  totalArea+=(wgt*detJ[id]);
	}
    }
}

void CutCellInterp(double *x, int d, int e, int p, double* Jinv, 
   		   double *ijk, double *xcut, 
		   double* rst) //, double *invJcut, double *detJcut)
//
// Routine that outputs bases and basis derivs at a specified rst 
// location on the cut cell. Used for both face and volume quad pts
//
// Inputs: 
//   x = global coords of orig tri
//   d = number of spat dims
//   e = original element type
//   Jinv = orig element inv Jac
//   ijk = cut cell rst value of desired pt
//   xcut = global coords of cut region
//
// Outputs
//   rst = rst value of desired point
{
 
  // notation: 
  //   xyz = global coord sys
  //   rst = local coord sys of original cell
  //   ijk = local coord sys of cut cell

  int b,w,i,j,l,ld,ij,m,n;
  int nbasisx=order2basis[e][1]; 
  int nbasis=order2basis[e][p]; 
  double jcut[d*d],mat[d][d],ijac[d][d],det;
  double xycut[d];
  double bd[nbasis][d]; 
  l=m=ld=ij=0;
  det = -1000; 


  //get xyz coord of nodes on the two elements
  //corresponding to rst = 0,0,1

  double x0[2] = {0.0, 0.0}; 
  double x0cut[2] = {0.0, 0.0}; 
  double u[2] = {0,0};
  for(i=0;i<nbasisx;i++){
    x0cut[0] = x0cut[0] + xcut[i]*basis[e][i](u);
    x0cut[1] = x0cut[1] + xcut[i+nbasisx]*basis[e][i](u);
  }
  for(i=0;i<nbasis;i++){
    x0[0] = x0[0] + x[i]*basis[e][i](u);
    x0[1] = x0[1] + x[i+nbasis]*basis[e][i](u);
  }

  // Compute p=1 jacobian for subcell element, dx/di
  // only p=1 required for this regardless of element p 
  JacP1Tri(jcut,xcut,&det); 

  // Compute the global x y location of this cut cell quad pt
  // [xycut] = [jac_cut]*[ij] + x0_cut
  axb(jcut,ijk,xycut,d);
  for(j=0;j<d;j++) xycut[j] = xycut[j] + x0cut[j] - x0[j];
  //printf("x0    : %f %f\n",x0[0],x0[1]);
  //printf("xycut : %f %f\n",xycut[0],xycut[1]);
  // Compute rst coord of cut cell quad
  // [rs] = [Jinv_orig][xycut - x0_orig]
  axb(Jinv,xycut,rst,d);

  if(rst[0]<-1e-14 || rst[1]<-1e-14 || rst[0] > 1+1e-10 || rst[1] > 1+1e-10){
    printf("ERROR: rst = %f %f\n",rst[0],rst[1]);

    printf("\tx0 = (%f %f),",x0[0],x0[1]); 
    u[0] = 1; 
    u[1] = 0;
    x0[0]=x0[1] = 0.0;  
    for(i=0;i<nbasis;i++){
      x0[0] = x0[0] + x[i]*basis[e][i](u);
      x0[1] = x0[1] + x[i+nbasis]*basis[e][i](u);
    }
    printf("\t(%f %f)",x0[0],x0[1]); 
    u[0] = 0; 
    u[1] = 1;
    x0[0]=x0[1] = 0.0;  
    for(i=0;i<nbasis;i++){
      x0[0] = x0[0] + x[i]*basis[e][i](u);
      x0[1] = x0[1] + x[i+nbasis]*basis[e][i](u);
    }
    printf("\t(%f %f)\n",x0[0],x0[1]); 


    printf("\tx0cut = (%f %f),",x0cut[0],x0cut[1]); 
    u[0] = 1; 
    u[1] = 0;
    x0cut[0]=x0cut[1] = 0.0;  
    for(i=0;i<nbasisx;i++){
      x0cut[0] = x0cut[0] + xcut[i]*basis[e][i](u);
      x0cut[1] = x0cut[1] + xcut[i+nbasisx]*basis[e][i](u);
    }
    printf("\t(%f %f),",x0cut[0],x0cut[1]); 
    u[0] = 1; 
    u[1] = 0;
    x0cut[0]=x0cut[1] = 0.0;  
    for(i=0;i<nbasisx;i++){
      x0cut[0] = x0cut[0] + xcut[i]*basis[e][i](u);
      x0cut[1] = x0cut[1] + xcut[i+nbasisx]*basis[e][i](u);
    }
    printf("\t(%f %f)\n",x0cut[0],x0cut[1]); 
    printf("\torig Jinv = %f %f ; %f %f\n",Jinv[0],Jinv[1],Jinv[2],Jinv[3]);
     
    exit(1);  
  }
}

void BasesVCut(double *x, double *Jinv,double *detJ,
               double *xcut, double *bvcut, double *bvdcut, 
               double *Jinvcut, double *detJcut,
               int d, int e, int p,int icut)

// Routine that builds up the bases and derivs of the interior cut cell quad pts
{
  int g = p2g[e][p];
  int l,n,m,ld,ii,w,b,i,j; 
  int nbasis1=order2basis[e][1];       
  int nbasis=order2basis[e][p+(p==0)];
  double u[2],bd[nbasis][d],mat[d][d],jac[d][d],det;

  //Loop over vol quad pts  
  double ijk[d]; 
  l=n=m=ld=ii=0;
  for(w=0; w<ngElem[e][p]; w++){

    // get local coord of quad pt
    for(j=0; j<d;j++) ijk[j] = gauss[e][g][(d+1)*w+j]; //double check inputs to gauss XXX

    // get bases and derivs at quad pt
    CutCellInterp(x,d,e,p,Jinv,ijk,xcut,u); 

    //store bvcut at each quad pt
    for(b=0;b<nbasis;b++){ // loop over bases
      for(j=0;j<d;j++){
	    bd[b][j]=basis_d[e][b*d+j](u); // accumulate bases derivs at quad pt
      }
      if (p > 0) bvcut[l++]=basis[e][b](u); // filled in as bv[nGL][nbasis]
    }
    if (p==0) bvcut[l++]=1;



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
//	   bvdcut[ld]=jac[0][i]*bd[b][0];
	   bvdcut[ld]=Jinv[i]*bd[b][0];
	   for(j=1;j<d;j++)
//	     bvdcut[ld]+=jac[j][i]*bd[b][j]; 
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

void CutFaceWeights(double *x, double *Jinv, int pc, int* iptr, double *xcut, double *bfcutL, double *bfdcutL, double *bfcutR, double *bfdcutR, double *Jinvcut, double *faceWeight, int d, int e, int p, int* cut2neigh, int iorig, int icut)
{
  int b,w,i,j,ij,l,ld,n,m,f,f1,eid,neigh;
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

  m=l=0;
  ld=ij=0;
  for(f=0;f<nfaces;f++){    
      // for every Gauss-point on this sub face
      for(w=0;w<ngGL[e][p];w++){

  	    v=gaussgl[e][g][(d)*w];  // gauss location from 0 to 1
	    wgt=gaussgl[e][g][(d)*w+1];	 // gauss weight

  	    // this is specific to 2-D 
	    // remember fill order will be reversed for element
            // sharing this edge   
            f1=(f+1)%nfaces;

	    // do 1D interp to get rst coord (on sub element) of quad pt
	    ijk[0]=(1-v)*eloc[e][p][d*f]  +v*eloc[e][p][d*f1]; 
	    ijk[1]=(1-v)*eloc[e][p][d*f+1]+v*eloc[e][p][d*f1+1];

            // now convert ijk coord on sub elem to rst coord on full L side elem
            eid = iorig; 
            CutCellInterp(x+iptr[pc*eid+1],d,e,p,Jinv+iptr[pc*eid+8],ijk,xcut,uL); //,Jinvcut,detJcut);  
	    
	    eid = cut2neigh[f];
            if(eid>-1){
              CutCellInterp(x+iptr[pc*eid+1],d,e,p,Jinv+iptr[pc*eid+8],ijk,xcut,uR); //,Jinvcut,detJcut);  
	    }

	    // get the bases and derivs at gauss legendre pts
	    for(b=0;b<nbasis;b++){
	      for(j=0;j<d;j++){
		bdL[b][j]=basis_d[e][b*d+j](uL); 
		if(eid!=-1)
		  bdR[b][j]=basis_d[e][b*d+j](uR); 
	      }
	      // basis[][](u) computes basis function at u
	      if (p > 0){
		bfcutL[l]=basis[e][b](uL);  // filled as bf[nfaces][nGL][nbasis]
		if(eid>-1){
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
  	    // XXX nbasis doesn't match xcut here
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

            // get [dx/dr]^-1
            if (d==2) invmat2x2(mat,jacL,det);

//Commenting all this out b/c we don't ever even use bd at faces
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
   } //loop over faces
}

void cut_mass_matrix(double *M, double *x, double *Jinv, double *xcut, double *detJcut, int d, int e, int p)
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
              CutCellInterp(x,d,e,p,Jinv,ijk,xcut,u); //,invJcut,detJcut); 
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
void mass_matrix(double *M, double *x, int d, int e, int p, int debug)
{
  int i,j,ij,w,b,ii,jj;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];
  double mat[d][d],jac[d][d];
  double u[d],wgt,det;
  int g=p2g[e][p+1]; // gauss point type for this element type, use one
                     // order higher making sure mass matrix is exact

  if (p > 0 ) {
    for(i=0;i<nbasis;i++)
      for(j=0;j<nbasis;j++)
	{
	  ij=nbasis*i+j;
	  M[ij]=0;
	  for(w=0;w<ngElem[e][p+1];w++) // 
	    {
	      for(jj=0;jj<d;jj++)
		u[jj]=gauss[e][g][(d+1)*w+jj];
	      /* evaluate the jacobian, it's not stored at this many locations (p=p+1)*/
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
	      wgt=gauss[e][g][(d+1)*w+2]*det;
	      /* could use bv here instead of reevaluating */
	      /* revaluate this if mesh is deforming */
	      M[ij]+=(wgt*basis[e][i](u)*basis[e][j](u));
	    }
	}
if(debug){
printf("\tfull Jac = %f %f %f %f\n",mat[0][0], mat[0][1], mat[1][0], mat[1][1]);
printf("\tfull invJac = %f %f %f %f\n",jac[0][0], jac[0][1], jac[1][0], jac[1][1]);
}	      
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

void Jacobian(double *x,double *bv, double *bvd, double *Jinv, 
              double *detJ, int d, int e, int p)
{
  int b,w,i,j,l,ld,ij,m,n;
  double u[d];
  double mat[d][d],jac[d][d],det;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];  // basis derivative  
  int g=p2g[e][p];  // gauss quadrature data for this element type
  //double xx,yy;
  l=n=m=ld=ij=0;
  for(w=0;w<ngElem[e][p];w++) // loop over gauss pts
    {
      for(j=0;j<d;j++) // get rs coordinates of gauss pt
	u[j]=gauss[e][g][(d+1)*w+j];
      //xx=yy=0;
      for(b=0;b<nbasis;b++) // loop over bases
	{
	  for(j=0;j<d;j++) 
	    bd[b][j]=basis_d[e][b*d+j](u); // accumulate bases derivs at quad pt
	  if (p > 0) bv[l++]=basis[e][b](u); // filled in as bv[nGL][nbasis]
	}
      if (p==0) bv[l++]=1;

      for(b=0;b<nbasis;b++)
        for(i=0;i<d;i++)

      for(b=0;b<nbasis;b++){
}
      //build jacobian dx/dr
      for(i=0;i<d;i++)
	{
	  for(j=0;j<d;j++)
	    {
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
	for(j=0;j<d;j++) {
	  Jinv[ij++]=jac[i][j];
        }

      // get basis derivs dN/dx = dN/dr*(dx/dr)^-1
      if (p > 0) {
      for(b=0;b<nbasis;b++)
	for(i=0;i<d;i++)
	  {
	     bvd[ld]=jac[0][i]*bd[b][0];
	      for(j=1;j<d;j++)
		bvd[ld]+=jac[j][i]*bd[b][j];
	     ld++;
	  }
      }
      else {
        for(i=0;i<d;i++) bvd[ld++]=0;
      }
      detJ[n++]=det;
    }
}

void FaceWeights(double *x, double *bf, double *bfd, double *Jinv, double *faceWeight, 
		 int d, int e, int p)
{
  int b,w,i,j,ij,l,ld,n,m,f,f1;
  double u[d],v,wgt;
  double mat[d][d],jac[d][d],det;
  int nbasis=order2basis[e][p+(p==0)];
  double bd[nbasis][d];  // basis derivative  
  int nfaces=facePerElem[e];
  int g=p2gf[e][p]; // gauss quadrature type for this element type
  double xx,yy;
  // ok, we are not going above 3D now
  double Ja[3];
  double Jb[3]={0,0,1};

  // for every face
  m=l=0;
  ld=ij=0;

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
	  // u is coordinate of gauss pt in rst
	  u[0]=(1-v)*eloc[e][p][d*f]  +v*eloc[e][p][d*f1]; //eloc is node location in rst; doing 1D interp?
	  u[1]=(1-v)*eloc[e][p][d*f+1]+v*eloc[e][p][d*f1+1];
	  //xx=yy=0;
	  
	  // get the bases and derivs at gauss legendre pts
	  for(b=0;b<nbasis;b++)
	    {
	      for(j=0;j<d;j++)
		{
		  bd[b][j]=basis_d[e][b*d+j](u); 
		  //bfd[ld++]=bd[b][j];      // filled in as bfd[nfaces][nGL][nbasis][d];
		}
	      // basis[][](u) computes basis function at u
	      if (p > 0) bf[l++]=basis[e][b](u);  // filled as bf[nfaces][nGL][nbasis]
	      //xx+=x[b]*bf[l-1];
	      //yy+=x[nbasis+b]*bf[l-1];
	    }
	  if (p==0) bf[l++]=1.0;
	  //printf("(r,s,xx,yy)=%f %f %f %f\n",r,s,xx,yy);

	  // build the jacobians of the edge
	  for(i=0;i<d;i++)
	    {
              faceWeight[d*m+i]=0;
	      Ja[i]=0;
	      for(j=0;j<d;j++)
		{
	          mat[i][j]=x[i*nbasis]*bd[0][j];
		  for(b=1;b<nbasis+(nbasis==1);b++)
		    mat[i][j]+=x[i*nbasis+b]*bd[b][j]; 
		} 
	      /* need a cross product here */
              for(j=0;j<d;j++)
	       {
                Ja[i]+=(mat[i][j]*face2elem[e][d*f+j]); // face2elem goes from edge coord to rst elem coord
		//TODO add Jb[i] calculation for 3D elements here
	       }
            }
	    // do faceWeight = Ja x zhat
	    // This gives me the normal vector
	    cross(&(faceWeight[2*m]),Ja,Jb,d); 
            //printf("fw:(%f %f)\n",faceWeight[2*m],faceWeight[2*m+1]);

            // get [dx/dr]^-1
            if (d==2) invmat2x2(mat,jac,det);

            // compute [dr/dx][dN/dr]
	    if (p > 0) {
 	     for(b=0;b<nbasis;b++)
	      {
		for(i=0;i<d;i++)
		  {
		    bfd[ld]=jac[0][i]*bd[b][0];
		    for(j=1;j<d;j++)
		       bfd[ld]+=jac[j][i]*bd[b][j];
		    ld++;
		  }
	      }
	    }
	    else {
              for(i=0;i<d;i++) bfd[ld++]=0;
	    }
	    for(i=0;i<d;i++)
              for(j=0;j<d;j++)
	          Jinv[ij++]=jac[i][j];
	    m++;
	}
    }
}

void COMPUTE_GRID_METRICS(double *x, double *bv, double *bvd,double *JinvV, 
			  double *detJ,double *bf, double *bfd, double *JinvF, double *faceWeight,
			  int *iptr, int d, int e, int p, int nelem, int pc, int* iblank, int imesh)
{
  int i,j,b;
  int ip,ix,ibv,ibvd,ibf,ij,idetj,ibfd,ijf,ifw;
  
  // Metrics for the uncut cells
  for(i=0;i<nelem;i++)
    {
      iblank[i] = 0; 

      ip=pc*i;
      ix   =iptr[ip+1];
      ibv  =iptr[ip+2];
      ibvd =iptr[ip+3];
      ij   =iptr[ip+4];
      idetj=iptr[ip+5];

      ibf  =iptr[ip+6];
      ibfd =iptr[ip+7];
      ijf  =iptr[ip+8];
      ifw  =iptr[ip+9];

//printf("ComputeGrid Metrics Elem %i\n", i); 
      Jacobian(x+ix, bv+ibv, bvd+ibvd, JinvV+ij,detJ+idetj,d,e,p); // basis on vol
      FaceWeights(x+ix,bf+ibf,bfd+ibfd,JinvF+ijf,faceWeight+ifw,d,e,p); // basis on face
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
			 int* cut2neigh, int imesh)
{
  int i,j,b,eid;
  int ip,ix,ij,idetj,ijf;
  int cip,cix,cibv,cibvd,cibf,cij,cidetj,cibfd,cijf,cifw,cc2n;

  // Metrics for the cut regions
  for(i=0;i<necut;i++)
  {
    eid = cut2e[i]; // Get elem ID of original elem 
    ip = eid*pc; // Get elem ID of original elem 

    // Quantities for the original element
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

    // get bases at cut vol quad pts
    BasesVCut(x+ix, JinvV+ij, detJ+idetj, 
              xcut+cix, bvcut+cibv, bvdcut+cibvd,
              JinvVcut+cij, detJcut+cidetj,
              d,e,p,i); 


    // basis on face
    CutFaceWeights(x, JinvF,pc,iptr,xcut+cix,
  		   bfcutL+cibf,bfdcutL+cibfd,
  		   bfcutR+cibf,bfdcutR+cibfd,
                   JinvFcut+cijf,fwcut+cifw,d,e,p,
		   cut2neigh+cc2n,eid,i); 
    printf("------\n");
  }

}

void MASS_MATRIX(double *mass,double *x, int *iptr, int d, int e, int p, int nelem, int pc, int imesh)
{
  int i;
  int ix,im,ixc;
  //compute mass matrix of all elements
  int nbasis=order2basis[e][p]; // XXX double check last argument, should be p=1 for cut cell
  int debug;
  for(i=0;i<nelem;i++)
    {
if(imesh==1 && i==225 & 1==0){
 debug = 1; 
printf("\n\nDEBUG: Mesh %i, cell %i:\n",imesh,i); 
}
else{ debug = 0;}

      ix=iptr[pc*i+1];
      im=iptr[pc*i+10];
      mass_matrix(&(mass[im]),&(x[ix]),d,e,p,debug); 

if(debug){
int m = 0; 
for(int j=0;j<nbasis;j++)
for(int k=0;k<nbasis;k++){
 printf("Mesh %i, cell %i, full mass[%i] = %f\n",imesh,i,m,mass[im+m]);
m++;
}
}

    }
}

void CUT_MASS_MATRIX(double *mass,double *x, double *Jinv, int *iptr, double *xcut, double *detJcut, int *iptrc, int d, int e, int p, int nelem, int pc, int pccut, int necut, int* cut2e, int imesh)
{
  int i,j,k,ld,eid;
  int ix,im,ij,ixc,id;

  int nbasis=order2basis[e][p]; // XXX double check last argument, should be p=1 for cut cell
  //subtract out cut portion
  for(i=0;i<necut;i++){
      eid = cut2e[i]; // find original element id 
      ix=iptr[pc*eid+1];
      ij=iptr[pc*eid+4];
      im=iptr[pc*eid+10];

      ixc=iptrc[pccut*i+1];
      id=iptrc[pccut*i+5]; 
      cut_mass_matrix(&(mass[im]),&(x[ix]),&(Jinv[ij]),&(xcut[ixc]),&(detJcut[id]),d,e,p);

if(imesh==1 && eid==225){
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

