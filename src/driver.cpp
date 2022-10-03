#include "dgsand.h"

int main(int argc, char *argv[])
{
  const int nmesh=argc-3; // # of input files = nmesh
  printf("NMESH = %i \n",nmesh); 

  double offset,x0;
  sscanf(argv[argc-2],"%lf",&offset); // last input is offset in grid units
  printf("offset=%f\n",offset);
  sscanf(argv[argc-1],"%lf",&x0); // last input is offset in grid units
  printf("x0=%f\n",x0);

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
//  double x0=8.1235;//1.75;
//  double x0=9.218750; // XXX this needs to be dynamically computed based on offset
//double x0=9.98750; // 
//  double x0=1.75;
  for(i=0;i<nmesh;i++) {
    sol[i].setup(argv[i+1],i,offset);// x0 will come out of here
    sol[i].init(i);
    sol[i].mass_matrix(i);
    sol[i].initTimeStepping(i);
  }

printf("\nx0 = %f, offset = %f\n",x0,offset);

  if(nmesh>1){
    for(i=0;i<nmesh;i++) {
printf("\n=================\nCUTTING MESH %i\n=================\n",i);
      sol[i].cut(x0,i);
      sol[i].cut_metrics(x0,i);
    }
    for(i=0;i<nmesh;i++) {
      B = 1-i; 
printf("\n=================\nENTERING OVERSET SETUP FOR MESH %i\n=================\n",i);      
      sol[i].setupOverset(sol[B].iptr,
		          sol[B].iptrc,
		   	  sol[B].x,
		   	  sol[B].xcut,
			  sol[B].JinvV,
			  sol[B].cut2e,
		          sol[B].necut,
			  sol[B].nelem);

    }
  }	  

  int nsteps=sol[0].getNsteps();
  int nsave=sol[0].getNsave();
  double dt=sol[0].getDt();
  const double rk[4]={0.25,8./15,5./12,3./4};
  for(i=0;i<nmesh;i++)
    sol[i].output(i,0);
 
  for(int n=1;n<=nsteps;n++) {
    // RK step 1
    for(i=0;i<nmesh;i++){
      // exchange overset flux information	    


      // Euler
/*      
      if(nmesh>1){
printf("==================\nEXCHANGING OVERSET FOR MESH %i Step %i, Euler \n===================\n",i,n);
        B = 1-i; 
        sol[i].exchangeOverset(sol[B].q, sol[B].iptr,i); 
      }
printf("==================\nCOMPUTING MESH %i Step %i, Euler \n===================\n",i,n);
      sol[i].computeRHS(sol[i].q,i);
    }
    for(i=0;i<nmesh;i++)
      {
        sol[i].update(sol[i].q,sol[i].q,dt);
      }
*/

      // RK
      if(nmesh>1){
        B = 1-i; 
        sol[i].exchangeOverset(sol[B].q, sol[B].iptr,i); 
      }
printf("==================\nCOMPUTING MESH %i Step %i, RK 1\n===================\n",i,n);
      sol[i].computeRHS(sol[i].q,i);
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
        sol[i].exchangeOverset(sol[B].qstar, sol[B].iptr, i); 
      }
printf("==================\nCOMPUTING MESH %i Step %i, RK 2 \n===================\n",i,n);
      sol[i].computeRHS(sol[i].qstar,i);
    }
    for(i=0;i<nmesh;i++)
      sol[i].update(sol[i].qstar,sol[i].q,rk[2]*dt);

    //RK step 3
    for(i=0;i<nmesh;i++){
      if(nmesh>1){
        B = 1-i; 
        sol[i].exchangeOverset(sol[B].q, sol[B].iptr, i); 
      }
printf("==================\nCOMPUTING MESH %i Step %i, RK 3\n===================\n",i,n);
      sol[i].computeRHS(sol[i].qstar,i);
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
  
    // Output data
    if (n%nsave==0) {
      for(i=0;i<nmesh;i++)
	sol[i].output(i,n);
    }
  }
}
