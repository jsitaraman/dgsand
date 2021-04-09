#define invmat2x2(mat,jac,det) { det=mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];     \
                                 jac[0][0]=mat[1][1]/det;jac[0][1]=-mat[0][1]/det; \
                                 jac[1][0]=-mat[1][0]/det;jac[1][1]=mat[0][0]/det; }

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
