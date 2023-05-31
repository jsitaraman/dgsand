#include<math.h>
#define invmat2x2(mat,jac,det) { det=mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];     \
                                 jac[0][0]=mat[1][1]/det;jac[0][1]=-mat[0][1]/det; \
                                 jac[1][0]=-mat[1][0]/det;jac[1][1]=mat[0][0]/det; }

double roundeps(double in,double eps) {
 return round(in/eps)*eps;
}

void matmult(double *A, double *B, double *C,int m,int n, int p){
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

void lu(double* a, double* l, double* u, int n)
{
  // modified from
  // https://www.tutorialspoint.com/cplusplus-program-to-perform-lu-decomposition-of-any-matrix

  int i = 0, j = 0, k = 0;
  for(int i = 0;i<n*n;i++){
    u[i] = 0.0;
    l[i] = 0.0;
  }

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (j < i)
         l[j*n+i] = 0;
         else {
            l[j*n+i] = a[j*n+i];
            for (k = 0; k < i; k++) {
               l[j*n+i] = l[j*n+i] - l[j*n+k] * u[k*n+i];
            }
         }
      }
      for (j = 0; j < n; j++) {
         if (j < i)
         u[i*n+j] = 0;
         else if (j == i)
         u[i*n+j] = 1;
         else {
            u[i*n+j] = a[i*n+j] / l[i*n+i];
            for (k = 0; k < i; k++) {
               u[i*n+j] = u[i*n+j] - ((l[i*n+k] * u[k*n+j]) / l[i*n+i]);
            }
         }
      }
   }
}

void fsub(double* L, double* y, double* b, int n){
  double sum; 
  for(int i=0;i<n;i++) y[i] = 0.0; 
//  for(int i=0;i<n*n;i++)  printf("\tL(%i) = %.16e;\n",i,L[i]);

  for(int i=0;i<n;i++){    
    y[i] = b[i];
    sum = 0; 
    for(int j=0;j<i;j++)
      sum += L[i*n+j]*y[j];
    y[i] = (b[i]-sum)/(L[i*n+i]);
  } 
}

void bsub(double* U, double* x, double* b, int n){
  double sum; 
  for(int i=0;i<n;i++) x[i] = 0.0; 
  
  for(int i=n-1;i>=0;i--){
    x[i] = b[i];

    sum = 0;
    for(int j=n;j>i;j--)
      sum +=U[i*n+j]*x[j]; 
    x[i] = (b[i]-sum)/(U[i*n+i]);
  }
}


// New solver to hopefully improve accuracy for sliver elements
void lusolve_reg(double* A, double* b, int n, int neq)
{
  double  x[n], y[n], L[n*n], U[n*n], btmp[2*n];

  // Fill regularized Areg = [A; lambda*I] and breg = [b; zeros(n,1)]
  double Areg[2*n*n], AregT[2*n*n], AregA[n*n], breg[2*n], Aregb[n];
  for(int i=0;i<2*n;i++){
    for(int j=0;j<n;j++){
      if(i>n){
        Areg[i*n+j] = 0.0;
      }
      else{
        Areg[i*n+j] = A[i*n+j];
      }
      AregT[j*n+i] = Areg[i*n+j];
    }
  }
  matmult(AregT,Areg,AregA,n,2*n,n);
  

  // Do LU solve
  lu(AregA,L,U,n);
/*
  for(int i=0;i<n;i++)
    if(L[i*n+i]==0 || U[i*n+i]==0){
      printf("\nFAIL: L or U has zero in diagonal!\n");
      for(int j=0;j<n;j++)
      for(int k=0;k<n;k++)
        printf("\nL(%i,%i) = %.16e;\n",j+1,k+1,L[j*n+k]);
      for(int j=0;j<n;j++)
      for(int k=0;k<n;k++)
        printf("\nU(%i,%i) = %.16e;\n",j+1,k+1,U[j*n+k]);
    }
*/
  for(int i=0;i<neq;i++){
    for(int j=0;j<2*n;j++)
      if(j<n){
        breg[j] = b[n*i+j];     
      } 
      else{
        breg[j] = 0.0;
      } 
    matmult(AregT,breg,Aregb,n,2*n,1);
    fsub(L,y,Aregb,n);
    bsub(U,btmp,y,n); 
    for(int j=0;j<n;j++) b[n*i+j] = btmp[j];
  }
}

void lusolve(double* A, double* b, int n, int neq)
{
  double  x[n], y[n], L[n*n], U[n*n], btmp[n];
  double eps = 1e-15; 

  lu(A,L,U,n);
/*
  for(int i=0;i<n;i++)
    if(L[i*n+i]==0 || U[i*n+i]==0){
      printf("\nFAIL: L or U has zero in diagonal!\n");
      for(int j=0;j<n;j++)
      for(int k=0;k<n;k++)
        printf("\nL(%i,%i) = %.16e;\n",j+1,k+1,L[j*n+k]);
      for(int j=0;j<n;j++)
      for(int k=0;k<n;k++)
        printf("\nU(%i,%i) = %.16e;\n",j+1,k+1,U[j*n+k]);
      exit(1); 
    }
  for(int i=0;i<n;i++)
  for(int j=0;j<n;j++)
    printf("L(%i,%i) = %.16e;\n",i+1,j+1,L[i*n+j]);
  for(int i=0;i<n;i++)
  for(int j=0;j<n;j++)
    printf("U(%i,%i) = %.16e;\n",i+1,j+1,U[i*n+j]);
  */

  for(int i=0;i<neq;i++){
//    for(int j=0;j<n;j++) btmp[j] = roundeps(b[n*i+j],eps);
    for(int j=0;j<n;j++) btmp[j] = b[n*i+j];
    fsub(L,y,btmp,n);
    bsub(U,btmp,y,n); // rewriting b (aka res) with final solution
    for(int j=0;j<n;j++) b[n*i+j] = btmp[j];
  }
}

void checksol(double* A, double* x, double* b, int n, int ind)
{
  // checks the accuracy of the computed solution and for the existence of NaNs

  // first check solution for 
  for(int i=0; i<n; i++)
    if(isnan(x[i])){
      printf("\nERROR: NaN found in solution of element %i. Exiting.\n",ind);
      exit(1); 
    }

  // get residual of linear solve by doing
  // res = Ax-b
  double *res; 
  res=(double *)calloc(n,sizeof(double));
  matmult(A,x,res,n,n,1);
  for(int i=0;i<n;i++) 
    res[i] = res[i]-b[i];

/*
  if(debug){
    for(int i=0;i<n;i++) 
//      printf("\tAx-b[%i] = %.16e\n",i,res[i]);

    if(fabs(res[0])>1e-13){

      for(int i = 0;i<n;i++)
      for(int j = 0;j<n;j++)
        printf("debugA(%i,%i) = %.16e;\n",i+1,j+1, A[i*n+j]);

      for(int i = 0;i<n;i++)
        printf("debugx(%i) = %.16e;\n",i+1, x[i]);

      for(int i = 0;i<n;i++)
        printf("debugb(%i) = %.16e;\n",i+1, b[i]);
    }
  }
*/

  free(res);  
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

  // Replace matrix with Areg'*Areg and vec with Areg'*b
  debug2=0;
  if(debug) debug2 = 0; 
  matmult(AregT,Areg,AregA,n,2*n,n);
  if(debug) debug2 = 1; 
  matmult(AregT,breg,Aregb,n,2*n,neq);

for(i=0;i<neq;i++)
  for(j=0;j<n;j++){
    ind1 = i*n+j;    
    ind2 = j*neq+i;
    b_in[ind1] = Aregb[ind2]; 
  }

  a=(double **)malloc(sizeof(double)*n);
  double (*b)[n]=(double(*)[n]) b_in ;

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
  free(b);
  free(Areg);
  free(breg);
  free(AregA);
  free(Aregb);
  free(AregT);

  return;

}
