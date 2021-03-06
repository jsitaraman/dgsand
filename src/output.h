void output_coords(FILE *fp,double *x,double *q, int d, int e , int p, int nfields)
{
  int i,j, k,b;
  double u[d], bi;
  double xv[d],qv[nfields];
  int nbasis=order2basis[e][p];
  
  for(i=0;i<nbasis;i++)
    {
      for(j=0;j<nbasis;j++)
	{
	  u[0]=eloc[e][p][i*d+0];
	  u[1]=eloc[e][p][i*d+1];
	  bi=basis[e][0](u);
	  for(k=0;k<d;k++)
	    xv[k]=bi*x[k*nbasis];
	  for(k=0;k<nfields;k++)
	    qv[k]=bi*q[k*nbasis];
	  for(b=1;b<nbasis;b++)
	    {
	      bi=basis[e][b](u);
	      for(k=0;k<d;k++)
		xv[k]+=(bi*x[k*nbasis+b]);
	      for(k=0;k<nfields;k++)
		qv[k]+=(bi*q[k*nbasis+b]);
	    }
	}
  for(k=0;k<d;k++)
    fprintf(fp,"%18.16f ",xv[k]);
  for(k=0;k<nfields;k++)
    fprintf(fp,"%18.16f ",qv[k]);
  fprintf(fp,"\n");
  }
}

void output_connectivity(FILE *fp,int offset,int p)
{
  int nsubtri[3]={1,1,4};
  int subtri_p1[3]={1,2,3};
  int subtri_p2[12]={1,4,6,4,2,5,6,4,5,6,5,3};
  int *subtri[3]={subtri_p1,subtri_p1,subtri_p2};
  int i;
  
  for(i=0;i<nsubtri[p];i++)
    fprintf(fp,"%d %d %d\n",subtri[p][3*i+0]+offset,subtri[p][3*i+1]+offset,subtri[p][3*i+2]+offset);
}

void OUTPUT_TECPLOT(int step,double *x, double *q,
		    int pc, int *iptr, int pde, int d, int e, int p, int nelem)
{

  char fname[80];
  char intstring[7];
  char qstr[3];
  int nfields=get_nfields[pde](d);
  int nvert=order2basis[e][p+(p==0)];
  int nsubtri[3]={1,1,4};
  
  FILE *fp;
  int ix,iq,i,j,k;
  //
  sprintf(intstring,"%d",100000+step);
  sprintf(fname,"flow%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"DGSAND output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\", ");
  for(i=0;i<nfields;i++)
    {
      sprintf(qstr,"Q%d",i);
      fprintf(fp,"\"%s\",",qstr);
    }
  fprintf(fp,"\n");
  if (p > 0) {
      fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=TRIANGLE, F=FEPOINT\n",nelem*nvert,nelem*nsubtri[p]);
  }
  else {
    fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=TRIANGLE, F=FEBLOCK\n",nelem*nvert,nelem*nsubtri[p]);
    fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL");
    for(i=0;i<nfields;i++) fprintf(fp, ", %d=CELLCENTERED",i+3);
    fprintf(fp,")\n");
  }
  if (p==0) {
    for (j=0;j<d;j++)
      for(i=0;i<nelem;i++)
	for(k=0;k<nvert;k++)
	  fprintf(fp,"%f\n",x[i*d*nvert+j*nvert+k]);
    for(j=0;j<nfields;j++)
      for(i=0;i<nelem;i++)
	fprintf(fp,"%f\n",q[i*nfields+j]);
  }
  else {
    for(i=0;i<nelem;i++)
      {
	ix=iptr[pc*i+1];
	iq=iptr[pc*i];
	output_coords(fp,x+ix,q+iq,d,e,p,nfields);
      }
  }
  for(i=0;i<nelem;i++)
    output_connectivity(fp,i*nvert,p);
  fclose(fp);
  
}
