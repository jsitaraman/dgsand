/* read a 2-D ascii grid of triangles of a given format
   make this adapt to various formats later, perhaps VTK */
void readgrid2D(double **xcoord, int **elem2node,int **ibc,
		int *p,  int *nnodes, int *nelem, int *nbnodes)
{
  int i,j,pg,m,total_n,n;
  int bcnode,bctype;
  char line[256];
  FILE *fp=fopen("grid.dat", "r");
  
  fgets(line,256,fp);  
  sscanf(line,"%d %d %d",nnodes,nelem,p);
  //printf("%d %d %d\n",*nnodes,*nelem,*p);
  (*xcoord)=(double *)malloc(sizeof(double)*(*nnodes)*2);
  for(i=0;i<(*nnodes);i++)
    {
      fgets(line,256,fp);
      sscanf(line,"%lf %lf",&((*xcoord)[2*i]),&((*xcoord)[2*i+1]));
    }

  pg=order2basis[0][*p+(*p==0)];
  (*elem2node)=(int *)malloc(sizeof(int)*(*nelem)*pg);

  for(i=0;i<(*nelem);i++)
  {
      fgets(line,256,fp);
      m=0;
      total_n=0;
      while (1 == sscanf(line + total_n, "%d%n", &j, &n))
	    {
	    total_n += n;
	    (*elem2node)[pg*i+m]=j-1;
	    m++;
	    if (m==pg) break;
	    }
  }
  fgets(line,256,fp);
  sscanf(line,"%d",nbnodes);
  (*ibc)=(int *)malloc(sizeof(int)*(*nbnodes)*2);
  for(i=0;i<(*nbnodes);i++)
   {
    fscanf(fp,"%d %d",&bcnode,&bctype);
    (*ibc)[2*i]=bcnode-1;
    (*ibc)[2*i+1]=bctype;
   }
  fclose(fp);
}
