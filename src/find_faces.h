// -*- c++ -*-
//
void insert_edge(int cellindx,int faceindx,
		 int e[2],int *flist, int *iptr, int *nfaces)
{
  int emin;
  int indx;
  int g[2];

  emin=e[0] < e[1] ? e[0] : e[1]; // returns the minimum node # of the two
  indx=iptr[emin]; 
  while(indx!=-1) // second time seeing this node?
    {
      g[0]=flist[7*indx]; // find orig edge first associated with node emin 
      g[1]=flist[7*indx+1];
      if (((g[0]==e[0]) && (g[1]==e[1])) || // if we're seeing this edge for the second time, 
	  ((g[1]==e[0]) && (g[0]==e[1])))   // store right side info, else update indx and hop out of loop
	{
	  // right side
	  flist[7*indx+4]=cellindx;
	  flist[7*indx+5]=faceindx;
	  return;
	}
      indx=flist[7*indx+6]; // 
    }
  flist[7*(*nfaces)]=e[0];
  flist[7*(*nfaces)+1]=e[1];
  flist[7*(*nfaces)+2]=cellindx; // left side
  flist[7*(*nfaces)+3]=faceindx;
  flist[7*(*nfaces)+4]=-1; // right side
  flist[7*(*nfaces)+5]=-1;
  flist[7*(*nfaces)+6]=iptr[emin];
  // set the pointer to >-1 to indicate we've already looked at it
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
	  int nv, // 3
		int nvmax) // nbasisx == 3?
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
  // loop through element vertices (3*nelem) and
  // accumulate the total number of unique nodes
  for(i=0;i<nv*nsurfcells;i++)
    // ternary operator = condition ? outcome1 : outcome2
    // bface = elem2node; elem2node gives node # based on elem #
    nnodes=nnodes > bface[i] ? nnodes : bface[i];
  iptr=(int *) malloc(sizeof(int)*nnodes);
  //
  *nfaces=0;
  // keeps track if it's the first time we touched this node
  for (i=0;i<nnodes;i++)
    iptr[i]=-1;
  
  bcnode=(int *)calloc(nnodes,sizeof(int));
  for(i=0;i<nbnodes;i++)
    bcnode[ibc[2*i]]=ibc[2*i+1];
  
  // loop over nelem
  for(i=0;i<nsurfcells;i++)
    {
      // loop over elem2nodes and add each edge
      for(n=0;n<nv;n++)
	{
	  np1=(n+1)%nv; // is just n+1 looping back to 0 
	  e[0]=bface[nvmax*i+n]; // node id 1, bface = elem2node
	  e[1]=bface[nvmax*i+np1]; // node id 2
	  insert_edge(i,n,e,flist,iptr,nfaces); // add edge to flist array
	}
    }
  // now we should have all of the edges and the cell/face ids for their L and R sides
  //
  // fill elem2face
  // map between all elems and their edges to the face id and its left and right sides 
  *elem2face=(int *) malloc(sizeof(int)*nsurfcells*nv); //3*nelem
  *faces=(int *) malloc(sizeof(int)*(*nfaces)*6); // 6 info per face, same as flist
  for(i=0;i<*nfaces;i++)
    {
      ileft=flist[7*i+2]; // cell index or -1 depending on orientation
      fleft=flist[7*i+3]; // face index or -1
      iright=flist[7*i+4]; // cell index or -1
      fright=flist[7*i+5]; // face index or -1 
      (*elem2face)[nv*ileft+fleft]=(i+1); //elem2face[3*cellID_l + faceID_l] = face id
      if (iright > -1) {
	      (*elem2face)[nv*iright+fright]=-(i+1);
      }
      else { // if this edge doesn't have a left and right side, it's a boundary node
        if (bcnode[flist[7*i]]==1 || bcnode[flist[7*i+1]]==1) flist[7*i+4]=-1;
        if (bcnode[flist[7*i]]==2 || bcnode[flist[7*i+1]]==2) flist[7*i+4]=-2;
        //printf("%d %d\n",i,flist[7*i+4]);
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

