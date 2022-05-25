#include<stdio.h>
double param[6];
void parseInputs(char *inputfile, char *gridfile,int *pde, int *itype, int *nsteps, double *dt, int *nsave, int *ireg)
{
  FILE *fp;
  char line[256];
  char comments[100];
  fp=fopen(inputfile,"r");
  fgets(line,256,fp);  sscanf(line,"pde=%d",pde);
  fgets(line,256,fp);  sscanf(line,"itype=%d",itype);
  fgets(line,256,fp);  sscanf(line,"nsteps=%d",nsteps);
  fgets(line,256,fp);  sscanf(line,"ireg=%d",ireg);
  fgets(line,256,fp);  sscanf(line,"dt=%lf",dt);
  fgets(line,256,fp);  sscanf(line,"nsave=%d",nsave);
  fgets(line,256,fp);  sscanf(line,"param[0]=%lf",&param[0]);
  fgets(line,256,fp);  sscanf(line,"param[1]=%lf",&param[1]);
  fgets(line,256,fp);  sscanf(line,"param[2]=%lf",&param[2]);
  fgets(line,256,fp);  sscanf(line,"param[3]=%lf",&param[3]);
  fgets(line,256,fp);  sscanf(line,"param[4]=%lf",&param[4]);
  fgets(line,256,fp);  sscanf(line,"param[5]=%lf",&param[5]);
  fgets(line,256,fp);  sscanf(line,"gridfile=%s",gridfile);
  fclose(fp);
  for(int i=0;i<6;i++) printf("%f ",param[i]);
}
void output_params()
{
  for(int i=0;i<6;i++) printf("%f ",param[i]);
}
/*
void main(void)
{
  int pde,itype,nsteps;
  double dt;
  parseInputs("input.dgsand",&pde,&itype,&nsteps,&dt);
  printf("%d %d %d %f\n",pde,itype,nsteps,dt);
  printf("%f %f %f %f %f %f\n",param[0],param[1],param[2],param[3],param[4],param[5]);
}
*/
