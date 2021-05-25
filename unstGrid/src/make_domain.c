/* 
program which generates input file for
the unstructured grid generator based on simple
domain description:
to be used in combination with ugrid.c which
is the unstructured grid generator.

Jayan Sitaraman
IIT Madras

last updated by J.Sitaraman
08/29/08

*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# define EPS 1E-04

/*

input file format

type = 1 (rectangle)
1 lx ly ux uy num direction
(lx,ly) bottom left, (ux,uy) upper right
num -> number of points
direction -> 1 for clockwise, 0 for anti-clockwise

type= 0 (circle)
0 cx cy R num direction
(cx,cy,R) -> center and radius, num and direction same as above

type = 2 (read from data file)
2 filename cx cy ang nsf direction
filename -> name of data file
cx,cy    -> offset to place the file
ang      -> angle of rotation
nsf -> user specified node spacing function, 0 for automatic calculation
direction same as previous cases

output file nodes

*/


main()
{
 FILE *fp,*fp1;
 char filename[40];
 int i,k=0,nb[10],numx,numy,num,num_domains,j,num_points;
 int type,direction;
 float psi,deltax,deltay,lx,ly,ux,uy,cx,cy,r,dtheta,theta,nsf1;
 float xval[5000],yval[5000],nsf[5000];
 float *xtemp,*ytemp;
 float xxc,yyc,ang;

 fp=fopen("domain","r");
 fscanf(fp,"%d",&num_domains);

 for(i=0;i<num_domains;i++)
   { 
     fscanf(fp,"%d",&type);
     if (type == 1)
       {
	 fscanf(fp,"%f %f %f %f %d %d",&lx,&ly,&ux,&uy,&num,&direction);

	 numx=num*(ux-lx)/(uy-ly+ux-lx)/2;
	 numy=(num-2*numx)/2;
	 if (direction==1)
	   {
	     deltax=2.0/(float)numx;
	     deltay=2.0/(float)numy;

	     for(psi=-1;psi<1-EPS;psi+=deltay)
	       {
		 yval[k]=(1-psi)*0.5*ly+(1+psi)*0.5*uy;
		 xval[k]=lx;
		 nsf[k]=0.5*deltay*(uy-ly);
		 k++;
	       }
	      for(psi=-1;psi<1-EPS;psi+=deltax)
	       {
		 xval[k]=(1-psi)*0.5*lx+(1+psi)*0.5*ux;
		 yval[k]=uy;
		 nsf[k]=0.5*deltax*(ux-lx);
		 k++;
	       }
	      for(psi=-1;psi<1-EPS;psi+=deltay)
	       {
		 yval[k]=(1-psi)*0.5*uy+(1+psi)*0.5*ly;
		 xval[k]=ux;
		 nsf[k]=0.5*deltay*(uy-ly);
		 k++;
	       }
	       for(psi=-1;psi<1-EPS;psi+=deltax)
	       {
		 xval[k]=(1-psi)*0.5*ux+(1+psi)*0.5*lx;
		 yval[k]=ly;
		 nsf[k]=0.5*deltax*(ux-lx);
		 k++;
	       }
	   }

	 if (direction==0)
	   {
	     deltax=2.0/(float)numx;
	     deltay=2.0/(float)numy;

	     for(psi=-1;psi<1-EPS;psi+=deltax)
	       {
		 xval[k]=(1-psi)*0.5*lx+(1+psi)*0.5*ux;
		 yval[k]=ly;
		 nsf[k]=0.5*deltax*(ux-lx);
		 k++;
	       }
	      for(psi=-1;psi<1-EPS;psi+=deltay)
	       {
		 yval[k]=(1-psi)*0.5*ly+(1+psi)*0.5*uy;
		 xval[k]=ux;
		 nsf[k]=0.5*deltay*(uy-ly);
		 k++;
	       }
	      for(psi=-1;psi<1-EPS;psi+=deltax)
	       {
		 xval[k]=(1-psi)*0.5*ux+(1+psi)*0.5*lx;
		 yval[k]=uy;
		 nsf[k]=0.5*deltax*(ux-lx);
		 k++;
	       }
	       for(psi=-1;psi<1-EPS;psi+=deltay)
	       {
		 yval[k]=(1-psi)*0.5*uy+(1+psi)*0.5*ly;
		 xval[k]=lx;
		 nsf[k]=0.5*deltay*(uy-ly);
		 k++;
	       }
	   }
	nb[i]=k;
       }
     if (type==0)
       {
	 fscanf(fp,"%f %f %f %d %d",&cx,&cy,&r,&num,&direction);
	 dtheta=2*M_PI/(float)num;
	 
	 if (direction==1)
	   {
	     for(theta=2*M_PI;theta>0;theta-=dtheta)
	       {
		 xval[k]=cx+r*cos(theta);
		 yval[k]=cy+r*sin(theta);
		 nsf[k]=r*dtheta;
		 k++;
	       }
	   }
	 if (direction==0)
	   {
	     for(theta=0;theta<2*M_PI;theta+=dtheta)
	       {
		 xval[k]=cx+r*cos(theta);
		 yval[k]=cy+r*sin(theta);
		 nsf[k]=r*dtheta;
		 k++;
	       }
	   }
	 nb[i]=k;
       }
     
     if (type==2)
       {
	 fscanf(fp,"%s %f %f %f %f %d",filename,&xxc,&yyc,&ang,
		&nsf1,&direction);
	 printf("%f %f %f %f\n",xxc,yyc,ang,nsf1);
	 ang*=M_PI/180.0;
	 fp1=fopen(filename,"r");
	 fscanf(fp1,"%d",&num_points);
	 xtemp=(float *)malloc(sizeof(float)*num_points);
	 ytemp=(float *)malloc(sizeof(float)*num_points);
	 for(j=0;j<num_points;j++)
	   fscanf(fp1,"%f %f",&xtemp[j],&ytemp[j]);
	 fclose(fp1);
	 if (direction==1)
	   {
	     for(j=0;j<num_points-1;j++)
	       {
		 xval[k]=xtemp[j]*cos(ang)-ytemp[j]*sin(ang)+xxc;
		 yval[k]=ytemp[j]*cos(ang)+xtemp[j]*sin(ang)+yyc;
		 nsf[k]=(nsf1!=0)?nsf1:sqrt((xtemp[j+1]-xtemp[j])
					    *(xtemp[j+1]-xtemp[j])+
					    (ytemp[j+1]-ytemp[j])
					    *(ytemp[j+1]-ytemp[j]))*0.5;
		 k++;
	       }
	     j=num_points-1;
	     xval[k]=xtemp[j]*cos(ang)-ytemp[j]*sin(ang)+xxc;
	     yval[k]=ytemp[j]*cos(ang)+xtemp[j]*sin(ang)+yyc;
	     nsf[k]=nsf[k-1];
	     k++;
	   }

	 if (direction==0)
	   {
	     for(j=num_points-1;j>0;j--)
	       {
		 xval[k]=xtemp[j]*cos(ang)-ytemp[j]*sin(ang)+xxc;
		 yval[k]=ytemp[j]*cos(ang)+xtemp[j]*sin(ang)+yyc;
		 nsf[k]=(nsf1!=0)?nsf1:sqrt((xtemp[j-1]-xtemp[j])
					    *(xtemp[j-1]-xtemp[j])+
					    (ytemp[j-1]-ytemp[j])
					    *(ytemp[j-1]-ytemp[j]));
		 k++;
	       }
	     j=0;
	     xval[k]=xtemp[j]*cos(ang)-ytemp[j]*sin(ang)+xxc;
	     yval[k]=ytemp[j]*cos(ang)+xtemp[j]*sin(ang)+yyc;
	     nsf[k]=nsf[k-1];
	     k++;
	   }
	 nb[i]=k;
       }
   }

fp1=fopen("nodes","w");
fprintf(fp1,"%d",k);
fprintf(fp1," %d\n",num_domains);
for(i=0;i<num_domains-1;i++)
fprintf(fp1," %d\n",nb[i]);
for(i=0;i<k;i++)
fprintf(fp1,"%f %f %f\n",xval[i],yval[i],nsf[i]);

}
  

