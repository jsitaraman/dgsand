#define gamma 1.4
#define nparam 20

double field_var[nparam];
double fsv[nparam];

int newtonian_fluid_get_nfields(int d)
{
  return (d==2) ? 4 : 5;
}
  

void newtonian_fluid_init(double *X, double *Q, int d, int nbasis, int itype)
{
  /* hard coded de, this needs to be read from an input file */
  double rinf=param[0];
  double uinf=param[1];
  double vinf=param[2];
  double pinf=param[3];
  double x0=param[4];
  double y0=param[5];

  double pi=4*atan(1.0);
  double cinf=sqrt(gamma*pinf/rinf);
  double sinf=pinf/pow(rinf,gamma);
  double tinf=pinf/rinf;
  double strnth=1.0;
  double sigma=1.0; // 0.020; //1.0;
  double scale=1.0;
  double gm1=gamma-1;
  double afac=strnth/(2.0*pi);
  double bfac = -0.5*gm1/gamma*afac*afac;
  double xx,yy,rsq,ee,u,v,tpr,t,p,rho;
  
  int i;

  if (d==2) 
    {
      fsv[0]=rinf;
      fsv[1]=uinf*rinf;
      fsv[2]=vinf*rinf;
      fsv[3]=pinf/(gamma-1)+0.5*(rinf*(uinf*uinf+vinf*vinf));

      if (itype==0) {
	for(i=0;i<nbasis;i++)
	  {
	    xx=X[2*i];
	    yy=X[2*i+1];

	    Q[0 +     i]=rinf; // +0.01*(xx*xx+yy*yy); //*exp(-(xx*xx+yy*yy));
	    Q[nbasis+ i]=Q[i]*uinf;
	    Q[2*nbasis + i]=Q[i]*vinf;
	    Q[3*nbasis + i]=pinf/(gamma-1)+(0.5*Q[i]*(uinf*uinf+vinf*vinf));
	  }
      }
      if (itype==1) {
	/* isentropic vortex */
	for(i=0;i<nbasis;i++)
	  {
	    xx=X[2*i]  - x0;
	    yy=X[2*i+1] - y0;
	    rsq=xx*xx+yy*yy;
	    ee=exp(0.5*(1-sigma*rsq))*scale;
	    u=uinf-(afac*ee)*yy;
	    v=vinf+(afac*ee)*xx;
	    tpr=bfac*ee*ee;
	    t=tinf+tpr;
	    p=pow(pow(t,gamma)/sinf,1./gm1);
	    rho=p/t;	      
	    Q[0 +     i]=rho; 
	    Q[nbasis+ i]=Q[i]*u;
	    Q[2*nbasis + i]=Q[i]*v;
	    Q[3*nbasis + i]=p/(gamma-1)+(0.5*Q[i]*(u*u+v*v));
	  } 
      }
    } 
}
 
void newtonian_fluid_flux2D(double *F, double *Q,double *QD, int idir)
{
  double vel=Q[idir+1]/Q[0];
  double faceSpeed=0.0; // for later when we have ALE terms
  double rvel=vel-faceSpeed;
  double pressure = (Q[3]-0.5*(Q[2]*Q[2]+Q[1]*Q[1])/Q[0])*(gamma-1);
  /* TODO use QD for viscous flux */

  F[0]=rvel*Q[0];
  F[1]=rvel*Q[1];
  F[2]=rvel*Q[2];
  F[3]=rvel*Q[3]+pressure*vel;
  F[idir+1]=F[idir+1]+pressure;
}

void newtonian_fluid_farfield_values(double *F)
{
  F[0]=fsv[0];
  F[1]=fsv[1];
  F[2]=fsv[2];
  F[3]=fsv[3];
}

void newtonian_fluid_inviscid_wall_bc(double *FR, double *FL,double *normal, int d)
{
  double unorm[d],u[d],udotn,umag;

  unorm[0]=normal[0];
  unorm[1]=normal[1];
  umag=sqrt(unorm[0]*unorm[0]+unorm[1]*unorm[1]);
  unorm[0]/=umag;
  unorm[1]/=umag;
  u[0]=FL[1]/FL[0];
  u[1]=FL[2]/FL[0];
  udotn=u[0]*unorm[0]+u[1]*unorm[1];

  FR[0]=FL[0];
  FR[1]=FL[0]*(u[0]-2*udotn*unorm[0]);
  FR[2]=FL[0]*(u[1]-2*udotn*unorm[1]);
  FR[3]=FL[3];
}

#include "roeflx.h"

typedef int (*GET_nfields)();
typedef void (*INIT_fields)();
typedef void (*FLUX_function)();
typedef void (*INTERFACE_flux)();
	      
GET_nfields get_nfields[1]={&newtonian_fluid_get_nfields};
INIT_fields init_fields[1]={&newtonian_fluid_init};
FLUX_function flux_function[1]={&newtonian_fluid_flux2D};
INTERFACE_flux gradient_indep_flux[1]={&roeflx};
INIT_fields far_field[1]={&newtonian_fluid_farfield_values};
INIT_fields  wall_bc[1]={&newtonian_fluid_inviscid_wall_bc};

int number_of_fields(int pde, int d) { return get_nfields[pde](d);}
