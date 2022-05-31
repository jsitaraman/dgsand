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
  else if(nmesh>2){
    printf("dgsand: Only a maximum of 2 grids are supported\n"); 
    exit(0);
  }
  dgsand *sol=new dgsand[nmesh];
  
  int i, B; 
  double x0=1.5;  
  for(i=0;i<nmesh;i++) {
    sol[i].setup(argv[i+1]);
    sol[i].init();
    sol[i].mass_matrix();
    sol[i].initTimeStepping(i);
  }

  if(nmesh>1){
    for(i=0;i<nmesh;i++) {
printf("\n=================\nCUTTING MESH %i\n=================\n",i);
      sol[i].cut(x0,i);
      sol[i].cut_metrics(x0);

printf("\n ENTERING OVERSET SETUP\n");      
      B = 1-i; 
      sol[i].setupOverset(sol[B].iptr,
		          sol[B].iptrc,
		   	  sol[B].x,
			  sol[B].JinvV,
		          sol[B].nelem);

    }
  }	  

  int nsteps=sol[0].getNsteps();
  int nsave=sol[0].getNsave();
  double dt=sol[0].getDt();
  const double rk[4]={0.25,8./15,5./12,3./4};
  
  for(int n=1;n<=nsteps;n++) {
    // RK step 1
    for(i=0;i<nmesh;i++){
      // exchange overset flux information	    
      if(nmesh>1){
        B = 1-i; 
printf("Exchanging for mesh %i\n",i);
        sol[i].exchangeOverset(sol[B].q, sol[B].iptr); 
      }
      sol[i].computeRHS(sol[i].q);
    }
    for(i=0;i<nmesh;i++)
      {
	sol[i].update(sol[i].qstar,sol[i].q,rk[1]*dt);
	sol[i].update(sol[i].q,sol[i].q,rk[0]*dt);
      }

    // RK step 2
    for(i=0;i<nmesh;i++){
      if(nmesh>1){
        B = 1-i; 
//       sol[i].exchangeOverset(sol[B].qstar, sol[B].iptr); 
      }
      sol[i].computeRHS(sol[i].qstar);
    }
    for(i=0;i<nmesh;i++)
      sol[i].update(sol[i].qstar,sol[i].q,rk[2]*dt);

    //RK step 3
    for(i=0;i<nmesh;i++){
      if(nmesh>1){
        B = 1-i; 
//       sol[i].exchangeOverset(sol[B].qstar, sol[B].iptr); 
      }
      sol[i].computeRHS(sol[i].qstar);
    }
    for(i=0;i<nmesh;i++)
      sol[i].update(sol[i].q,sol[i].q,rk[3]*dt);

    // compute norms
    int imax;
    double rmax,rnorm;
    for(i=0;i<nmesh;i++)
      {
	sol[i].rnorm(imax,rmax,rnorm,rk[3]*dt);
	printf("mesh%d : step %d\t%18.16f\t%d\t%18.16f\n",i,n,rnorm,imax,rmax);
      }
    if (n%nsave==0) {
      for(i=0;i<nmesh;i++)
	sol[i].output(i,n);
    }
  }
}
