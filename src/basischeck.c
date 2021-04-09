#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "solvec.h"
#include "basislib.h"
#include "quadrature.h"
double basis2D(double r, double s,int elemtype, int k)
{ 
 double phi;
 switch(elemtype)  {
  case 0:
     phi=1;
  case 1:
     if (k==0) phi=1-r-s;
     if (k==1) phi=r;
     if (k==2) phi=s;
  case 2:
     if (k==0) phi=1-r-s;
     if (k==1) phi=r;
     if (k==2) phi=s;
     if (k==3) phi=r*(1-r-s);
     if (k==4) phi=r*s;
     if (k==5) phi=s*(1-r-s);
 }
 return phi;
}


void main(void)
{
  int iflag;
  double r[6]={0,1,0,0.5,0.5,0.0};
  double s[6]={0,0,1,0.0,0.5,0.5};
  double x[6]={0,1,0.5,0.5,0.75,0.25};
  double y[6]={0,0,sqrt(3)/2,0.,sqrt(3)/4,sqrt(3)/4};
  double **M;
  double **N;
  double xx,yy,zz,rr,ss;
  int i,j;

  M=(double **)malloc(sizeof(double *)*6);
  N=(double **)malloc(sizeof(double *)*3);

  for(i=0;i<6;i++) M[i]=(double *)malloc(sizeof(double)*6);
  for(i=0;i<3;i++) N[i]=(double *)malloc(sizeof(double)*6);

  for(i=0;i<6;i++)
    {
      for(j=0;j<6;j++) {
	/*printf("%f ",basis2D(r[i],s[i],2,j));*/
	//M[i][j]=basis2D(r[i],s[i],2,j);
	M[i][j]=tri_basis[j](r[i],s[i]);
	printf("%f ",M[i][j]);
      }
      N[0][i]=x[i];
      N[1][i]=y[i];
      N[2][i]=x[i]*x[i]+2*x[i]*y[i]+y[i]*y[i];
      printf("\n");
    }
  printf("\n");
  solvec(M,N,&iflag,6,3);
  for(i=0;i<6;i++)
    printf("%f %f %f\n",N[0][i],N[1][i],N[2][i]);
  printf("\n");
  
  rr=1./3;
  ss=1./3;
  xx=yy=zz=0;
  for(i=0;i<6;i++) {
   xx+=(basis2D(rr,ss,2,i)*N[0][i]);
   yy+=(basis2D(rr,ss,2,i)*N[1][i]);
   zz+=(basis2D(rr,ss,2,i)*N[2][i]);
  }
  printf("%f %f %f %f\n",xx,yy,xx*xx+2*xx*yy+yy*yy,zz);
  printf("\n");
  for(int g=0;g<3;g++)
    {
      printf("{ ");
      for(j=0;j<ngGL[g];j++)
	printf(" %.15e, %.15e,\n",(1+gaussgl[g][2*j])*0.5,(gaussgl[g][2*j+1])*0.5);
      printf("}\n");
    }

}
