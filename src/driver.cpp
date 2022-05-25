#include "dgsand.h"

int main(int argc, char *argv[])
{
  const int nmesh=argc-1;
  if (nmesh == 0) {
   printf("dgsand: Need at least one input file as argument\n");
   printf("e.g. for a two mesh case\n");
   printf("$../src/dgsand input.dgsand1 input.dgsand2 \n");
   exit(0);
  }
  dgsand *sol=new dgsand[nmesh];
  
  double x0=0.5;  
  for(int i=0;i<nmesh;i++) {
    sol[i].setup(argv[i+1]);
    sol[i].init();
    //sol[i].cut(x0);
    sol[i].mass_matrix();
    //sol[i].cut_metrics(x0);
    sol[i].initTimeStepping();
  }

  int nsteps=sol[0].getNsteps();
  int nsave=sol[0].getNsave();
  double dt=sol[0].getDt();
  const double rk[4]={0.25,8./15,5./12,3./4};
  
  for(int n=1;n<=nsteps;n++) {
    // RK step 1
    for(int i=0;i<nmesh;i++)      
      sol[i].computeRHS(sol[i].q);
    for(int i=0;i<nmesh;i++)
      {
	sol[i].update(sol[i].qstar,sol[i].q,rk[1]*dt);
	sol[i].update(sol[i].q,sol[i].q,rk[0]*dt);
      }
    // RK step 2
    for(int i=0;i<nmesh;i++)
      sol[i].computeRHS(sol[i].qstar);
    for(int i=0;i<nmesh;i++)
      sol[i].update(sol[i].qstar,sol[i].q,rk[2]*dt);
    //RK step 3
    for(int i=0;i<nmesh;i++)
      sol[i].computeRHS(sol[i].qstar);
    for(int i=0;i<nmesh;i++)
      sol[i].update(sol[i].q,sol[i].q,rk[3]*dt);
    // compute norms
    int imax;
    double rmax,rnorm;
    for(int i=0;i<nmesh;i++)
      {
	sol[i].rnorm(imax,rmax,rnorm,rk[3]*dt);
	printf("mesh%d : step %d\t%18.16f\t%d\t%18.16f\n",i,n,rnorm,imax,rmax);
      }
    if (n%nsave==0) {
      for(int i=0;i<nmesh;i++)
	sol[i].output(i,n);
    }
  }
}
