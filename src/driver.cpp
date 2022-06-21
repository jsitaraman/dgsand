#include "dgsand.h"

int main(int argc, char *argv[])
{
  const int nmesh=argc-2;
  printf("NMESH = %i \n",nmesh); 
  double ds;
  sscanf(argv[argc-1],"%lf",&ds);
  printf("ds=%f\n",ds);
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
  //double x0=9.9; //68750000000000-0*.625;
  double x0=10-ds*0.1; //-0*.625;
  for(i=0;i<nmesh;i++) {
    sol[i].setup(argv[i+1]);
    /*
    if (i==1) {
      sol[i].transform(1.8,0);
    }
    */
    sol[i].init(i);
    sol[i].mass_matrix(i);
  }

  if(nmesh>1){
    for(i=0;i<nmesh;i++) {
//printf("\n=================\nCUTTING MESH %i\n=================\n",i);
      sol[i].cut(x0,i);
      sol[i].cut_metrics(x0,i);

//printf("\n ENTERING OVERSET SETUP\n");      
      B = 1-i; 
      sol[i].setupOverset(sol[B].iptr,
		          sol[B].iptrc,
		   	  sol[B].x,
			  sol[B].JinvV,
		          sol[B].nelem);

    }
  }	  
  for(i=0;i<nmesh;i++)
    sol[i].initTimeStepping(i);

  int nsteps=sol[0].getNsteps();
  int nsave=sol[0].getNsave();
  double dt=sol[0].getDt();
  const double rk[4]={0.25,8./15,5./12,3./4};
  for(i=0;i<nmesh;i++)
    sol[i].output(i,0);
  double cons=0.0;
  for(i=0;i<nmesh;i++)
   cons+=sol[i].cons_metric(0);
  printf("cons : %.16e\n",cons);
  double cons0=cons;
  int euler=0;
  for(int n=1;n<=nsteps;n++) {
    // RK step 1
    for(i=0;i<nmesh;i++){
      // exchange overset flux information	          
      // Euler
      if(nmesh>1){
        B = 1-i; 
        sol[i].exchangeOverset(sol[B].q, sol[B].iptr,i); 
      }
      sol[i].computeRHS(sol[i].q,i);
    }
    if (euler) {
      for(i=0;i<nmesh;i++)
	{
	  sol[i].update(sol[i].q,sol[i].q,dt);
	}
    }
    else {
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
	//printf("==================\nCOMPUTING MESH %i Step %i, RK 2 \n===================\n",i,n);
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
	//printf("==================\nCOMPUTING MESH %i Step %i, RK 3\n===================\n",i,n);
	sol[i].computeRHS(sol[i].qstar,i);
      }
      for(i=0;i<nmesh;i++)
	sol[i].update(sol[i].q,sol[i].q,rk[3]*dt);
    }

    // compute norms
    int imax;
    double rmax,rnorm;
    cons=0;
    for(i=0;i<nmesh;i++)
      {
	sol[i].rnorm(imax,rmax,rnorm,rk[3]*dt);
	cons+=sol[i].cons_metric(0);
	printf("mesh%d : step %d\t%18.16f\t%d\t%18.16f\n",i,n,rnorm,imax,rmax);
      }
    printf("cons : %.16e %.16e %.16e\n",cons0,cons,abs(cons-cons0));
    if (n%nsave==0) {
      for(i=0;i<nmesh;i++)
	sol[i].output(i,n);
    }
    printf("l2 error : ");
    for(i=0;i<nmesh;i++)
      printf("%0.16e ",sol[i].compute_error(i,n));
    printf("\n");
  }
}
