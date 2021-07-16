// -*- c++ -*-
//
void insert_edge(int cellindx,int faceindx,
		 int e[2],int *flist, int *iptr, int *nfaces)
{
  int emin;
  int indx;
  int g[2];

  emin=e[0] < e[1] ? e[0] : e[1];
  indx=iptr[emin];
  while(indx!=-1)
    {
      g[0]=flist[7*indx];
      g[1]=flist[7*indx+1];
      if (((g[0]==e[0]) && (g[1]==e[1])) ||
	  ((g[1]==e[0]) && (g[0]==e[1])))
	{
	  flist[7*indx+4]=cellindx;
	  flist[7*indx+5]=faceindx;
	  return;
	}
      indx=flist[7*indx+6];
    }
  flist[7*(*nfaces)]=e[0];
  flist[7*(*nfaces)+1]=e[1];
  flist[7*(*nfaces)+2]=cellindx;
  flist[7*(*nfaces)+3]=faceindx;
  flist[7*(*nfaces)+4]=-1;
  flist[7*(*nfaces)+5]=-1;
  flist[7*(*nfaces)+6]=iptr[emin];
  iptr[emin]=*nfaces;
  (*nfaces)++;
}
//
void find_faces(int *bface,
		int **elem2face,
		int **faces,
		int *nfaces,
    int *ibc,
		int nsurfcells,
    int nbnodes,
	  int nv,
		int nvmax)
{
  int i,n,np1;
  int e[2];
  int *flist;
  int *iptr,*bcnode;
  int nnodes;
  int ileft,iright;
  int fleft,fright;
  //
  flist=(int *) malloc(sizeof(int)*nsurfcells*nv*7);
  nnodes=-1;
  for(i=0;i<nv*nsurfcells;i++)
    nnodes=nnodes > bface[i] ? nnodes : bface[i];
  iptr=(int *) malloc(sizeof(int)*nnodes);
  //
  *nfaces=0;
  //
  for (i=0;i<nnodes;i++)
    iptr[i]=-1;
  
  bcnode=(int *)calloc(nnodes,sizeof(int));
  for(i=0;i<nbnodes;i++)
    bcnode[ibc[2*i]]=ibc[2*i+1];
  
  //
  for(i=0;i<nsurfcells;i++)
    {
      for(n=0;n<nv;n++)
	{
	  np1=(n+1)%nv;
	  e[0]=bface[nvmax*i+n];
	  e[1]=bface[nvmax*i+np1];
	  insert_edge(i,n,e,flist,iptr,nfaces);
	}
    }
  //
  *elem2face=(int *) malloc(sizeof(int)*nsurfcells*nv);
  *faces=(int *) malloc(sizeof(int)*(*nfaces)*6);
  //
  for(i=0;i<*nfaces;i++)
    {
      ileft=flist[7*i+2];
      fleft=flist[7*i+3];
      iright=flist[7*i+4];
      fright=flist[7*i+5];
      (*elem2face)[nv*ileft+fleft]=(i+1);
      if (iright > -1) {
	      (*elem2face)[nv*iright+fright]=-(i+1);
      }
      else {
        if (bcnode[flist[7*i]]==1 || bcnode[flist[7*i+1]]==1) flist[7*i+4]=-1;
        if (bcnode[flist[7*i]]==2 || bcnode[flist[7*i+1]]==2) flist[7*i+4]=-2;
        printf("%d %d\n",i,flist[7*i+4]);
      }
      for(n=0;n<6;n++) (*faces)[6*i+n]=flist[7*i+n];
    }
  //
  //trace(*nfaces);
  //
  free(iptr);
  free(flist);
  free(bcnode);
}

