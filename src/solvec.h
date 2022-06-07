#define invmat2x2(mat,jac,det) { det=mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];     \
                                 jac[0][0]=mat[1][1]/det;jac[0][1]=-mat[0][1]/det; \
                                 jac[1][0]=-mat[1][0]/det;jac[1][1]=mat[0][0]/det; }

void matmult(double *A, double *B, double *C,int m,int n, int p, int debug){
// computes C = A*B where A(m,n) and B(n,p) and C(m,p)

  int i, j, k;
  int indA, indB, indC; 

  for(i=0;i<m;i++)
    for(j=0;j<p;j++){
      indC = i*p+j;
      C[indC] = 0.0; 
      for(k=0;k<n;k++){
        indA = i*n+k;
        indB = k*p+j;
        C[indC] += A[indA]*B[indB]; 
      }
    }

/*
 if(debug){
  for(i=0;i<m;i++)
    for(j=0;j<n;j++){
      indC = i*n+j;       
      printf("A_in(%i,%i) = %.8e;\n",i+1,j+1,A[indC]);
    }
  for(i=0;i<n;i++)
    for(j=0;j<p;j++){
      indC = i*p+j;       
      printf("B_in(%i,%i) = %.8e;\n",i+1,j+1,B[indC]);
    }
  for(i=0;i<m;i++)
    for(j=0;j<p;j++){
      indC = i*p+j; 
      printf("C_out(%i,%i) = %.8e;\n",i+1,j+1,C[indC]);
    }
//  exit(1); 
}
*/

}

void solvec(double **a,double **b,int *iflag,int n,int neq)
{
  int i,j,k,l,flag,temp1,m;
  double fact;
  double temp;
  double sum;
  double eps=1e-8;

  
  for(i=0;i<n;i++)
    {
      if (fabs(a[i][i]) < eps)
	{
	  flag=1;
	  for(k=i+1;k<n && flag;k++)
	    {
	      if (a[k][i]!=0)
                {
		  flag=0;
		  for(l=0;l<n;l++)
		    {
		      temp=a[k][l];
		      a[k][l]=a[i][l];
		      a[i][l]=temp;
		    }
		  for(m=0;m<neq;m++)
		    {
		      temp=b[m][k];
		      b[m][k]=b[m][i];
		      b[m][i]=temp;
		    }
                }
	    }
	  if (flag) {*iflag=0;return;}
	}
      for(k=i+1;k<n;k++)
	{
	  if (i!=k)
	    {
	      fact=-a[k][i]/a[i][i];
	      for(j=0;j<n;j++)
		{
		  a[k][j]+=fact*a[i][j];
		}
	      for(m=0;m<neq;m++)
		b[m][k]+=fact*b[m][i];
	    }
	}
    }

  for(i=n-1;i>=0;i--)
    {
      for(m=0;m<neq;m++)
	{
	  sum=0;
	  for(j=i+1;j<n;j++)
	    sum+=a[i][j]*b[m][j];
	  b[m][i]=(b[m][i]-sum)/a[i][i];
	}
    }
  *iflag=1;
  return;

}

void solvec_copy_reshape(double *a_in,double *b_in,int *iflag,int n,int neq)
{
  int i,j,k,l,flag,temp1,m;
  double fact;
  double temp;
  double sum;
  double eps=1e-8;
  double **a;
  double (*b)[n]=(double(*)[n])b_in;
  
  a=(double **)malloc(sizeof(double)*n);
  for(i=0;i<n;i++)
    {
      a[i]=(double *)malloc(sizeof(double)*n);
      for(j=0;j<n;j++)
	a[i][j]=a_in[i*n+j];
    }
  
  for(i=0;i<n;i++)
    {
      if (fabs(a[i][i]) < eps)
	{
	  flag=1;
	  for(k=i+1;k<n && flag;k++)
	    {
	      if (a[k][i]!=0)
                {
		  flag=0;
		  for(l=0;l<n;l++)
		    {
		      temp=a[k][l];
		      a[k][l]=a[i][l];
		      a[i][l]=temp;
		    }
		  for(m=0;m<neq;m++)
		    {
		      temp=b[m][k];
		      b[m][k]=b[m][i];
		      b[m][i]=temp;
		    }
                }
	    }
	  if (flag) {*iflag=0;return;}
	}
      for(k=i+1;k<n;k++)
	{
	  if (i!=k)
	    {
	      fact=-a[k][i]/a[i][i];
	      for(j=0;j<n;j++)
		{
		  a[k][j]+=fact*a[i][j];
		}
	      for(m=0;m<neq;m++)
		b[m][k]+=fact*b[m][i];
	    }
	}
    }

  for(i=n-1;i>=0;i--)
    {
      for(m=0;m<neq;m++)
	{
	  sum=0;
	  for(j=i+1;j<n;j++)
	    sum+=a[i][j]*b[m][j];
	  b[m][i]=(b[m][i]-sum)/a[i][i];
	}
    }
  *iflag=1;
  
  for(i=0;i<n;i++) free(a[i]);
  free(a);
  return;

}
void solvec_copy_reshape_reg(double *a_in,double *b_in,int *iflag,int n,int neq, int debug)
{
  int i,j,k,l,flag,temp1,m;
  double fact;
  double temp;
  double sum;
  double eps=1e-8;
  double **a;
  double lambda = 0.00; // controls how much regularization to do (lambda = 0 is smallest) 
  int ind1,ind2,debug2;
 
  // Fill regularized Areg = [A; lambda*I] and breg = [b; zeros(n,1)]
  double *Areg, *breg, *AregA, *Aregb, *AregT;
  Areg=(double *)calloc(2*n*n,sizeof(double));
  AregT=(double *)calloc(2*n*n,sizeof(double));
  breg=(double *)calloc(2*n*neq,sizeof(double)); // handle neq dim correctly
  AregA=(double *)calloc(n*n,sizeof(double));
  Aregb=(double *)calloc(n*neq,sizeof(double)); // handle neq dim correctly

/*
if(debug){
for(i=0;i<n;i++)
for(j=0;j<n;j++){
      ind1 = i*n+j;
printf("A_in(%i,%i) = %f\n",i+1,j+1,a_in[ind1]);
}
printf("\n");
for(j=0;j<n;j++)
printf("b_in(%i) = %f\n",j+1,b_in[j]);
printf("\n");
}
*/

  for(i=0;i<2*n;i++){
    for(j=0;j<n;j++){
      ind1 = i*n+j;
      if(i<n){
        Areg[ind1] = a_in[ind1];
      }
      else if((i-n)==j){
        Areg[ind1] = lambda; 
      }
    }
  
    if(i<n){
      for(k=0;k<neq;k++) breg[i*neq+k] = b_in[k*n+i];
///breg[i*neq+k] = b_in[k*n+i]; 
//breg[k*2*n+i] = b_in[k*n+i]; 
    }
    else{
      for(k=0;k<neq;k++) breg[i*neq+k] = 0.0; 
    }
  }
  for(i = 0; i<2*n; i++)
    for(j = 0; j<n; j++){
      ind1 = i*n+j; 
      ind2 = j*2*n+i; 
      AregT[ind2] = Areg[ind1]; 
    }

/*
if(debug){
for(i=0;i<2*n;i++)
for(j=0;j<n;j++){
ind1 = i*n+j;
printf("A_reg(%i,%i) = %f;\n",i+1,j+1,Areg[ind1]);
}
printf("\n");
for(j=0;j<2*n;j++)
for(k=0;k<neq;k++){
ind1 = j*neq+k;
printf("b_reg(%i,%i) = %f;\n",j+1,k+1,breg[ind1]);
}
printf("\n");
for(i=0;i<n;i++)
for(j=0;j<2*n;j++){
ind1 = i*2*n+j;
printf("A_regT(%i,%i) = %f;\n",i+1,j+1,AregT[ind1]);
}
}
*/
  // Replace matrix with Areg'*Areg and vec with Areg'*b
  debug2=0;
  if(debug) debug2 = 0; 
  matmult(AregT,Areg,AregA,n,2*n,n,debug2);
  if(debug) debug2 = 1; 
  matmult(AregT,breg,Aregb,n,2*n,neq,debug2);

for(i=0;i<neq;i++)
  for(j=0;j<n;j++){
    ind1 = i*n+j;    
    ind2 = j*neq+i;
    b_in[ind1] = Aregb[ind2]; 
  }

/*
 * if(debug){
for(i=0;i<n;i++)
for(j=0;j<n;j++){
ind1 = n*i+j; 
printf("AregA(%i,%i) = %f;\n",i+1,j+1,AregA[ind1]);
}
printf("\n");
for(j=0;j<n;j++)
printf("Aregb(%i) = %f;\n",j+1,Aregb[j]);
printf("\n");
}
*/

  a=(double **)malloc(sizeof(double)*n);
  double (*b)[n]=(double(*)[n]) b_in ;

/*
if(debug){
for(j=0;j<neq;j++){
for(i=0;i<n;i++){
printf("b_rsz(%i,%i) = %f;\n",j+1,i+1,b[j][i]); 
}}
}
*/
  for(i=0;i<n;i++)
    {
      a[i]=(double *)malloc(sizeof(double)*n);
      for(j=0;j<n;j++)
	a[i][j]=AregA[i*n+j];
    }

  for(i=0;i<n;i++)
    {
      if (fabs(a[i][i]) < eps)
	{
	  flag=1;
	  for(k=i+1;k<n && flag;k++)
	    {
	      if (a[k][i]!=0)
                {
		  flag=0;
		  for(l=0;l<n;l++)
		    {
		      temp=a[k][l];
		      a[k][l]=a[i][l];
		      a[i][l]=temp;
		    }
		  for(m=0;m<neq;m++)
		    {
		      temp=b[m][k];
		      b[m][k]=b[m][i];
		      b[m][i]=temp;
		    }
                }
	    }
	  if (flag) {*iflag=0;return;}
	}
      for(k=i+1;k<n;k++)
	{
	  if (i!=k)
	    {
	      fact=-a[k][i]/a[i][i];
	      for(j=0;j<n;j++)
		{
		  a[k][j]+=fact*a[i][j];
		}
	      for(m=0;m<neq;m++)
		b[m][k]+=fact*b[m][i];
	    }
	}
    }

  for(i=n-1;i>=0;i--)
    {
      for(m=0;m<neq;m++)
	{
	  sum=0;
	  for(j=i+1;j<n;j++)
	    sum+=a[i][j]*b[m][j];
	  b[m][i]=(b[m][i]-sum)/a[i][i];
	}
    }
/*
if(debug){
for(j=0;j<n;j++)
printf("b_out(%i) = %f\n",j+1,b[j]);
printf("\n");
}
*/
  *iflag=1;
  
  for(i=0;i<n;i++) free(a[i]);
  free(a);
  free(Areg);
  free(breg);
  free(AregA);
  free(Aregb);
  free(AregT);

  return;

}
