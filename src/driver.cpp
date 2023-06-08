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
 
  // ============= 
  // SETUP MESH 
  // ============= 
  int i, B; 
  for(i=0;i<nmesh;i++) {
    sol[i].setup(argv[i+1],i,offset);// x0 will come out of here
    sol[i].init(i);
  }
  printf("\nx0 = %f, offset = %f\n",x0,offset);

  if(nmesh>1){
    for(i=0;i<nmesh;i++) {
      printf("\n=================\nCUTTING MESH %i\n=================\n",i);
      // setup conservative overset cut cells
      sol[i].cut(x0,i);		// find cut cells
      sol[i].cut_metrics(x0,i); // find cut bases and etc
     
      // cell merging routines
      sol[i].findCellMerge(i); 	 // check to see if any cells need merging
      sol[i].findParents(i);	 // find parents for merged cells
//      sol[i].cellagg_metrics(i); // recompute bases for merged cells
    }

    // setup overset gauss pts and bases
    for(i=0;i<nmesh;i++){
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

  // Compute mass matrix (including merged cells)
  for(i=0;i<nmesh;i++) {
    sol[i].mass_matrix(i);
    sol[i].initTimeStepping(i);
  }
	  
  int nsteps=sol[0].getNsteps();
  int nsave=sol[0].getNsave();
  double dt=sol[0].getDt();
  const double rk[4]={0.25,8./15,5./12,3./4};
  for(i=0;i<nmesh;i++)
    sol[i].output(i,0);

  // Compute initial conservation metrics
  double cons=0.0,tmp=0;
  for(i=0;i<nmesh;i++){
   tmp=sol[i].cons_metric(0);
   printf("tmp = %f\n");
   cons+=tmp;
  }
  printf("cons : %.16e\n",cons);
  double cons0=cons;
  
  // ============= 
  // RUN TIMESTEPS 
  // ============= 
  int euler = 1;
  for(int n=1;n<=nsteps;n++) {
    // Euler
    if(euler){
      for(i=0;i<nmesh;i++){
      
        if(nmesh>1){
        printf("================================================\n");
        printf("EXCHANGING OVERSET FOR MESH %i Step %i, Euler \n",i,n); 
        printf("================================================\n");
          B = 1-i; 
          sol[i].exchangeOverset(sol[B].q, sol[B].iptr,i); 
        }
	printf("==================\nCOMPUTING MESH %i Step %i, Euler \n===================\n",i,n);
        sol[i].computeRHS(sol[i].q,i);
      }
      for(i=0;i<nmesh;i++)      
        sol[i].update(sol[i].q,sol[i].q,dt);
    } // end of euler
    else{  // RK
      for(i=0;i<nmesh;i++){
        // exchange overset flux information	    

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
      for(i=0;i<nmesh;i++) sol[i].update(sol[i].qstar,sol[i].q,rk[2]*dt);

      //RK step 3
      for(i=0;i<nmesh;i++){
        if(nmesh>1){
          B = 1-i; 
          sol[i].exchangeOverset(sol[B].q, sol[B].iptr, i); 
        }
	printf("==================\nCOMPUTING MESH %i Step %i, RK 3\n===================\n",i,n);
        sol[i].computeRHS(sol[i].qstar,i);
      }
      for(i=0;i<nmesh;i++) sol[i].update(sol[i].q,sol[i].q,rk[3]*dt);
    } // euler or rk
    
    // compute norms
    printf("\nCOMPUTING NORMS\n"); 
    int imax;
    double rmax,rnorm;
    cons=0;
    for(i=0;i<nmesh;i++)
      {
        sol[i].rnorm(imax,rmax,rnorm,rk[3]*dt);
        cons+=sol[i].cons_metric(0);
        printf("mesh%d : step %d\t%18.16f\t%d\t%18.16f\n",i,n,rnorm,imax,rmax);
      }
      printf("cons : %.16e %.16e %.16e\n",cons0,cons,fabs(cons-cons0));
  
    // Output data
    if (n%nsave==0) 
      for(i=0;i<nmesh;i++)
        sol[i].output(i,n);

    printf("l2 error : ");
    for(i=0;i<nmesh;i++)
      printf("%0.16e ",sol[i].compute_error(i,n));
    printf("\n");

  } // timestep
}
