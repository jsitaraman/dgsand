/* 2-D unstructured grid generation code using
   Constraint Delaunay triangulation technique

   Algorithm is hybrid between that used by 
   Sloan S.W and Raveendra V.V.S (see references). But there are
   some subtle differences to improve efficiency and
   quality of triangulation.
  
   Author : Jayanarayanan Sitaraman
   Final project :: AS 668 (CFD class)
   started :: August 1996

   Last updated by j.sitaraman on
   08/29/08 (packaged with examples)
*/

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# define swp(x,y)  {x=x+y;y=x-y;x=x-y;}
# define EPS 1e-02
# define float double

/* stack for pushing nodes */
typedef struct STACK
{
  struct NODE *node;
  struct STACK *next;
} STACK;

/* node data structure containing coordinates
   index, node spacing function */
typedef struct NODE
{
  float x;
  float y;
  float x1;  
  float y1;
  int num;
  float nsf;
  struct NODE *next;
  struct STACK *nearnode;
} NODE;

/* triangle data structure containing coordinates
   index, and pointers of nearest triangles */

typedef struct TRIANGLE
{
  struct NODE *node[3];
  struct TRIANGLE *neighbour[3];
  struct TRIANGLE *next;
}  TRIANGLE;

/* edge data structure containing coordinates 
   index and pointers to the next edge in list */

typedef struct EDGE
{
  struct NODE *f;
  struct NODE *s;
  float x2;  
  float y2;
  struct EDGE *next;
} EDGE;

/* triangle stack */

typedef struct STACK_TR
{
  struct TRIANGLE *tri;
  struct STACK_TR *next;
} STACK_TR;

/* delaunay flip stack */

typedef struct DELSTACK
{
  TRIANGLE *tri,*tri1;
  struct DELSTACK *next;
}DELSTACK;


void getvalues();
void initialise ( NODE *node_head, EDGE *edge_head);
int checkedge(NODE *N1,NODE *N2,EDGE *edge_head);
int maketriangle(NODE *node_head, EDGE *edge_head ,
		 TRIANGLE *triangle_head);
int checkedge(NODE *N1,NODE *N2,EDGE *edge_head);
int check_intersection(NODE *N1,NODE *N2,EDGE *edge_head);
int left (NODE *N1,NODE *N2,NODE *N3);
float findlength(NODE *N1,NODE *N2,NODE *N3);
void print_values(int n, TRIANGLE *triangle_head);
void push_proper(STACK_TR *stackhead,TRIANGLE *triangle);
TRIANGLE * check_common(STACK_TR *stack1, STACK_TR *stack2,TRIANGLE *triangle);
void set_neighbours(TRIANGLE *triangle_head,int Nodenum);
int check_circle (NODE * n[4]);
int getcircumcircle(NODE *N1,NODE *N2,NODE *N3,float *cx,float *cy,float *r);
int make_delaunay1(TRIANGLE *triangle_head,int num_triangles);
void swap(TRIANGLE **a,TRIANGLE **b);
float distance(NODE *N1,NODE *N2);
int make_newtriangles(int num_triangles,TRIANGLE *triangle_head,
		      NODE *node_head);
int check_nsf(TRIANGLE *triangle);
void push(STACK *stack,NODE *N);
void delete(STACK *stack,NODE *N);
int smooth2(NODE *node_head);
void print_fast(int n, int number,TRIANGLE *triangle_head,NODE *node_head);
void make_junk(int n, TRIANGLE *triangle_head,NODE *node_head);
void insert_node(TRIANGLE *triangle,TRIANGLE *triangle_head,NODE *node);
void pop_tr(DELSTACK *stack_head,TRIANGLE **tri,TRIANGLE **tri1);
void push_tr(DELSTACK *stack_head,TRIANGLE *tri,TRIANGLE *tri1);
void make_delaunay(DELSTACK *stack_head);

float *xval,*yval,*nsf;
int *pin;
int *np;
int nnodesf,nelemf;
float *xout;
int *tri;
int number,NB;
int numedges;
int plot_flag=0;

void generate_grid_(float *XVAL,float *YVAL,float *NSF,int *NPDOM,int *NN,int *NDOM)
{
 EDGE *edge_head,*edge;
 NODE *node_head,*node;
 TRIANGLE *triangle_head,*triangle;
 STACK *stack_head,*stack;

 FILE *fp,*fp1;
 int i,n,j,m,temp,gx,gy,complete,iter,iref=0,smoothing;

 NB=number=*NN;
 printf("NB=%d\n",NB);
 xval=(float *)malloc(sizeof(float)*NB);
 yval=(float *)malloc(sizeof(float)*NB);
 nsf=(float *)malloc(sizeof(float)*NB);
 pin=(int *)malloc(sizeof(float)*(NB+1));
 np=(int *)malloc(sizeof(int)*(*NDOM));

 for(i=0;i<NB+1;i++) pin[i]=0;
 for(i=1;i<(*NDOM);i++)
  {
    pin[NPDOM[i-1]]=1;
    np[i-1]=NPDOM[i-1];
  }

 for(i=0;i<NB;i++)
  {
   xval[i]=XVAL[i];
   yval[i]=YVAL[i];
   nsf[i]=NSF[i];
  }

 edge_head=(EDGE *)malloc(sizeof(EDGE));
 node_head=(NODE *)malloc(sizeof(NODE));
 edge_head->next=NULL;
 node_head->next=NULL;
 triangle_head=(TRIANGLE *)malloc(sizeof(TRIANGLE));

 initialise( node_head , edge_head);
 n=maketriangle(node_head,edge_head,triangle_head);
 set_neighbours(triangle_head,number);
 
 iter=0;
 complete=1;
 while(complete && iter<100)
   { 
     complete=make_delaunay1(triangle_head,n);
     iter++;
     if (iter%10==0) printf("boundary delaunay iteration:%d\n",iter);
   }
       
 iter=0;
 complete=1;
 iref++;  
 temp=n;
 n=make_newtriangles1(n,triangle_head,node_head);


 smoothing=0;
 while(complete && iter<5)
   {
     if (smoothing==0)
       smooth2(node_head);
     else
       smooth3(node_head);
     iter++;
     complete=make_delaunay1(triangle_head,n);
   }
 
 printf("number of triangles= %d\n",n);
 printf("number of nodes    =%d\n",number);

 nnodesf=number;
 nelemf=n;
 xout=(float *)malloc(sizeof(float)*nnodesf*2);
 tri=(int *)malloc(sizeof(int)*nelemf*3);
 
 node=node_head;
 m=0;
 while(node->next!=NULL)
   {
     node=node->next;
     xout[m++]=node->x;
     xout[m++]=node->y;
   }

 triangle=triangle_head;
 m=0;
 while(triangle->next!=NULL)
   {
     triangle=triangle->next;
     if (left(triangle->node[0],triangle->node[1],triangle->node[2]))
	{
	  tri[m++]=triangle->node[0]->num+1;
          tri[m++]=triangle->node[1]->num+1;
          tri[m++]=triangle->node[2]->num+1;
        }
     else
	{
          tri[m++]=triangle->node[0]->num+1;
          tri[m++]=triangle->node[2]->num+1;
          tri[m++]=triangle->node[1]->num+1;
        }
    }

 /* free all memory */

 node=node_head;
 while(node_head!=NULL)
 {
  node_head=node->next;
  /*
  stack_head=node->nearnode->next;
  stack=stack_head;
  printf("%p %p %p\n",node,stack_head,stack->next);
  while(stack_head!=NULL)
    {
      printf("%p\n",stack->next);
      stack_head=stack->next;
      free(stack);
      stack=stack_head;
    }      
  */
  free(node);
  node=node_head;
 }

 edge=edge_head;
 while(edge_head!=NULL)
 {
  edge_head=edge->next;
  free(edge);
  edge=edge_head;
 }

 triangle=triangle_head;
 while(triangle_head!=NULL)
 {
  triangle_head=triangle->next;
  free(triangle);
  triangle=triangle_head;
 }

 free(xval);
 free(yval);
 free(nsf);
 free(pin); 
 free(np);

}

void get_elem_count_(int *NNODES, int *NELEM)
{
 *NNODES=nnodesf;
 *NELEM=nelemf;
}

void get_tess_(float *X, int *TRI)
{
 int i;
 for (i=0;i<2*nnodesf;i++) X[i]=xout[i];
 for (i=0;i<3*nelemf;i++) TRI[i]=tri[i];
 free(xout);
 free(tri);
 number=0;
}

void initialise ( NODE *node_head, EDGE *edge_head)
{
  int i,k=0;
  NODE *node,*tnum;
  EDGE *edge;

  node=node_head;
  edge=edge_head;
  
  for(i=0;i<number;i++)
    { 
      node->next=(NODE *)malloc(sizeof(NODE));
      node=node->next;
      node->next=NULL;
      node->num=i;
      node->x=xval[i];
      node->y=yval[i];
      node->nsf=nsf[i];
      node->nearnode=NULL;
    }

  node=node_head;
  tnum=node->next;
  for(i=0;i<number;i++)
    {
      node=node->next;
      if (pin[i+1]==0 && i<number-1)
	{
	  edge->next=(EDGE *)malloc(sizeof(EDGE));
	  edge=edge->next;
	  edge->next=NULL;
	  edge->f=node;
	  edge->s=node->next;
	  k++;
	}
      else
	{
	  edge->next=(EDGE*)malloc(sizeof(EDGE));
	  edge=edge->next;
	  edge->next=NULL; 
	  edge->f=node;
	  edge->s=tnum;
	  tnum=node->next;
	  k++;
	}
    }
 numedges=k;
}    


void getvalues()
{
 FILE *fp,*fp1;
 int i,num_boundaries,nb,val=0;
 

fp = fopen("nodes","r");
fscanf(fp,"%d" ,&number); 
NB=number;
fscanf(fp,"%d",&num_boundaries);

for(i=1;i<num_boundaries;i++)
  {
    fscanf(fp,"%d",&nb);
    pin[nb]=1;
    np[i-1]=nb;
  }


for(i=0;i<number;i++)
fscanf(fp,"%lf %lf %lf",&xval[i],&yval[i],&nsf[i]);

fclose(fp);
}

  
int maketriangle(NODE *node_head, EDGE *edge_head ,TRIANGLE *triangle_head)
{
  int numnod=numedges;
  int i;
  NODE *N1,*N2,*N3,*Ntemp;
  float temp,temp1;
  int val,flag,k=0,check_val;
  NODE *node;
  EDGE *edge,*tempedge;
  TRIANGLE *triangle;

  node=node_head;
  edge=edge_head;
  triangle=triangle_head;

  while(edge_head->next!=NULL )
   {
     edge=edge_head->next;
     N1=edge->f;
     N2=edge->s;
     flag =0;
     check_val=0;
     tempedge=edge_head;
     while(tempedge->next !=NULL)
       {
	 tempedge=tempedge->next;
	 Ntemp=tempedge->f;
	 if (left(N1,N2,Ntemp ))
	   if (check_intersection(N1,Ntemp,edge_head))
	     if (check_intersection(N2,Ntemp,edge_head))
	       {
		 check_val=1;
		 temp=findlength(N1,N2,Ntemp);
		 if (flag==0)
		   {
		     flag=1;
		     temp1=temp;
		     N3=Ntemp;
		   }
		 else
		   {
		     if (temp<temp1)
		       {
			 N3=Ntemp;
			 temp1=temp;
		       }
		   }
	       }
       }
     if (check_val==0)
       {
	 edge=edge_head->next;
	 printf("*********************************************\n");
	 printf("not found for edge %d-%d\n",edge->f->num,edge->s->num);
	 printf("*********************************************\n");
	 printf("sorry..I can't handle very close points, if you are using\n");
	 printf("airfoil data please check the leading edge points \n");
	 exit(0);
       }

     numnod--;
     tempedge=edge_head->next;
     edge_head->next=tempedge->next;
          
     if (checkedge(N2,N3,edge_head))
       {
	 tempedge=edge_head->next;
	 edge_head->next=(EDGE *)malloc(sizeof(EDGE));
	 if (left(N2,N3,N1))
	   {
	     edge_head->next->f=N3;
	     edge_head->next->s=N2;
	     edge_head->next->next=tempedge;
	   }
	 else
	   {
	     edge_head->next->f=N2;
	     edge_head->next->s=N3;
	     edge_head->next->next=tempedge;
	   }
	 numnod++;
       }
     else
       numnod--;
          
     if (checkedge(N1,N3,edge_head))
       {
	 tempedge=edge_head->next;
	 edge_head->next=(EDGE *)malloc(sizeof(EDGE));
	 if (left(N1,N3,N2))	 	 
	 {
	     edge_head->next->f=N3;
	     edge_head->next->s=N1;
	     edge_head->next->next=tempedge;
	   }
	 else
	   {
	     edge_head->next->f=N1;
	     edge_head->next->s=N3;
	     edge_head->next->next=tempedge;
	   }
	 numnod++;
       }
     else
       numnod--;
     
     triangle->next=(TRIANGLE *)malloc(sizeof(TRIANGLE));
     triangle=triangle->next;
     triangle->next=NULL;
     triangle->node[0]=N1;
     triangle->node[1]=N2;
     triangle->node[2]=N3;
     /*   printf("%d: %d %d %d\n",k,N1->num,N2->num,N3->num);*/
     k++;
   }
  return k;
}


int checkedge(NODE *N1,NODE *N2,EDGE *edge_head)
{
 int i,val=1;
 EDGE *edge,*prevedge;
 edge=edge_head;
 while(edge->next!=NULL && val)
   {
     prevedge=edge;
     edge=edge->next;
     if ((edge->f==N1 && edge->s==N2)  || (edge->s==N1 && edge->f==N2)) 
       {
	 val=0;
	 prevedge->next=edge->next;
       }
   }
 return val;
}



int check_intersection(NODE *N1,NODE *N2,EDGE *edge_head)
{
  float x,y,x1,y1,a,a1,c2,c1,denom;
  int i;
  EDGE *edge;

  x=N2->x-N1->x;
  y=N2->y-N1->y;

  edge=edge_head;
  while(edge->next!=NULL)
   {
     edge=edge->next;
     if (!(N1==edge->f || N1==edge->s || N2==edge->f
	 || N2==edge->s))
       {
	 x1=edge->s->x-edge->f->x;
	 y1=edge->s->y-edge->f->y;
	 
	 c1=edge->f->x-N1->x;
	 c2=edge->f->y-N1->y;
	 
	 denom= -x*y1+x1*y;
	 
	 if (denom!=0)
	   {
	     a=(x1*c2-y1*c1)/denom;
	     a1=-(y*c1-c2*x)/denom;
	     
	     if ( a>EPS && a<1-EPS && a1>EPS && a1<1-EPS)
	       {
		 return 0;
	       }
	   }
       }
   }
  return 1;
}


int left (NODE *N1,NODE *N2,NODE *N3)
{
 float x,y,x1,y1;


 x=N1->x-N3->x;
 y=N1->y-N3->y;

 x1=N2->x-N3->x;
 y1=N2->y-N3->y;

 if ((x==0 && y==0) || (x1==0 && y1==0))  return 0;


return (x*y1-y*x1 >0)?1:0;

}


float findlength(NODE *N1,NODE *N2,NODE *N3)
{
  float x,y;
  float x1,y1;

  x=N1->x-N3->x;
  y=N1->y-N3->y;

  x1=N2->x-N3->x;
  y1=N2->y-N3->y;

  return sqrt(x*x+y*y)+sqrt(x1*x1+y1*y1);
}

void print_values(int n, TRIANGLE *triangle_head)
{
 FILE *fp;
 TRIANGLE *triangle;
 int i,j;

 fp = fopen("check.dat","w");

 fprintf(fp,"%d\n",n);
 triangle=triangle_head;

 while(triangle->next!=NULL)
 {
   triangle=triangle->next;
   fprintf(fp,"%f %f %f\n",triangle->node[0]->x, triangle->node[0]->y,0.0);
   fprintf(fp,"%f %f %f\n",triangle->node[1]->x, triangle->node[1]->y,0.0);
   fprintf(fp,"%f %f %f\n",triangle->node[2]->x, triangle->node[2]->y,0.0);
 }
 fclose(fp);
}

void make_junk(int n, TRIANGLE *triangle_head,NODE *node_head)
{
 TRIANGLE *triangle;
 NODE *node;
 FILE *fp;

 fp=fopen("junk.plt","w");

 fprintf(fp,"TITLE =\"triangle file\"\n");
 fprintf(fp,"VARIABLES =\"X\", \"Y\", \"\n");
 fprintf(fp,"ZONE T=\"DELAUNAY\",N= %d ,E= %d ,ET=TRIANGLE F=FEPOINT\n",number,n);

 node=node_head;
 while(node->next!=NULL)
   {
     node=node->next;
     fprintf(fp,"%f %f\n",node->x,node->y);
   }

 triangle=triangle_head;
 while(triangle->next!=NULL)
   {
     triangle=triangle->next;
     if (left(triangle->node[0],triangle->node[1],triangle->node[2]))
       fprintf(fp,"%d %d %d\n",triangle->node[0]->num+1,triangle->node[1]->num+1,
	     triangle->node[2]->num+1);
     else
       fprintf(fp,"%d %d %d\n",triangle->node[0]->num+1,triangle->node[2]->num+1,
	       triangle->node[1]->num+1);
   }

 fclose(fp);
/* 
 fp=fopen("neighbours","w");
 triangle=triangle_head;

 while(triangle->next!=NULL)
   {
     triangle=triangle->next;
     fprintf(fp,"%d %d %d\n",triangle->neighbour[0],triangle->neighbour[1],
	     triangle->neighbour[2]);
   }
 fclose(fp);
*/
 
}
   


void set_neighbours(TRIANGLE *triangle_head,int Nodenum)
{
 int i,dummy;
 STACK_TR *stack,*tmp;
 TRIANGLE *triangle;
 NODE *node;
 
 stack = (STACK_TR *)malloc(sizeof(STACK_TR)*(Nodenum));

 /* A stack is initialised for each node */
 for(i=0;i<Nodenum;i++)
   stack[i].next=NULL;
      
 /** the particular triangle number is pushed into the stack of each node 
     push_proper means usual kind of push as in it doesn't check whether
     its already there or not
     */
 
 triangle=triangle_head;
 while(triangle->next!=NULL)
   {
     triangle=triangle->next;
     push_proper(&stack[triangle->node[0]->num],triangle);
     push_proper(&stack[triangle->node[1]->num],triangle);
     push_proper(&stack[triangle->node[2]->num],triangle);
   }

 /* this loop essentially takes stacks corresponding to adjacent nodes on 
    each triangle and finds the common entry, this value is updated to 
    neighbouring triangle index
    */

 triangle=triangle_head;
 while(triangle->next!=NULL)
  {
    triangle=triangle->next;
    triangle->neighbour[0]=check_common(&stack[triangle->node[0]->num],&stack[triangle->node[1]->num],triangle);
    triangle->neighbour[1]=check_common(&stack[triangle->node[1]->num],&stack[triangle->node[2]->num],triangle);
    triangle->neighbour[2]=check_common(&stack[triangle->node[2]->num],&stack[triangle->node[0]->num],triangle);
  }
}

void push_proper(STACK_TR *stackhead,TRIANGLE *triangle)
{
  STACK_TR *temp;

  temp=stackhead->next;
  stackhead->next=(STACK_TR *)malloc(sizeof(STACK_TR));
  stackhead->next->tri=triangle;
  stackhead->next->next=temp;
}

TRIANGLE * check_common(STACK_TR *stack1, STACK_TR *stack2,TRIANGLE *triangle)
{
  STACK_TR *s1;
  STACK_TR *s2;

  s1=stack1;

  while(s1->next!=NULL)
    {
      s1=s1->next;
      s2=stack2;
      while(s2->next!=NULL)
	{
	  s2=s2->next;
	  if (s1->tri==s2->tri && s1->tri !=triangle) return s1->tri;
	}
    }
  return NULL;
}

int make_delaunay1(TRIANGLE *triangle_head,int num_triangles)
{
  int i,j,k,found,num;
  NODE *n[4];
  TRIANGLE *n13,*n30,*n20,*n12,*index;
  int l[3],p[3],complete=0,swapped;
  int *pin1;
  TRIANGLE *triangle;
  STACK *stck;
  float nx,ny,nsf;

  pin1=(int *)malloc(sizeof(int)*num_triangles);

  for(i=0;i<num_triangles;pin1[i++]=0);

  triangle=triangle_head;
  while(triangle->next!=NULL)
    {
      triangle=triangle->next;
      l[0]=distance(triangle->node[0],triangle->node[1]);
      l[1]=distance(triangle->node[1],triangle->node[2]);
      l[2]=distance(triangle->node[2],triangle->node[0]);
      
      p[0]=0;p[1]=1;p[2]=2;
      
      for(j=0;j<3;j++)
	for(k=j+1;k<3;k++)
	  if (l[j]<l[k]) swp(p[j],p[k]);
      
      for(j=0;j<3;j++)
	{
	  swapped=0;
	  if (triangle->neighbour[p[j]]!=NULL && !swapped)
	    {
	      if (p[j]==0)
		{
		  n[2]=triangle->node[2];
		  n[0]=triangle->node[0];
		  n[1]=triangle->node[1];
		  
		  n12=triangle->neighbour[1];
		  n20=triangle->neighbour[2];
		}
	      if (p[j]==1)
		{
		  n[2]=triangle->node[0];
		  n[0]=triangle->node[1];
		  n[1]=triangle->node[2];
		  
		  n12=triangle->neighbour[2];
		  n20=triangle->neighbour[0];
		}
	      if (p[j]==2)
		{
		  n[2]=triangle->node[1];
		  n[0]=triangle->node[2];
		  n[1]=triangle->node[0];
		  
		  n12=triangle->neighbour[0];
		  n20=triangle->neighbour[1];
		}
	      
	      
	      index=triangle->neighbour[p[j]];
	      
	      if (index->neighbour[0]==triangle)
		{
		  n[3]=index->node[2];
		  
		  n13=index->neighbour[1];
		  n30=index->neighbour[2];
		      
		  if (index->node[0]==n[1])
		    swap(&n13,&n30);
		}
	      
	      if (index->neighbour[1]==triangle)
		{
		  n[3]=index->node[0];
		  
		  n13=index->neighbour[2];
		  n30=index->neighbour[0];
		  
		  if (index->node[1]==n[1])
		    swap(&n13,&n30);
		}
	      
	      if (index->neighbour[2]==triangle)
		{
		  n[3]=index->node[1];
		  
		  n13=index->neighbour[0];
		  n30=index->neighbour[1];
		  
		  if (index->node[2]==n[1])
		    swap(&n13,&n30);
		}
		  		  	     	      
	      if (check_circle(n))
		{
		  swapped=1;
		  complete=1;
		  
		  if (n[2]->num >NB-1)  push(n[2]->nearnode,n[3]);
		  if (n[3]->num >NB-1)  push(n[3]->nearnode,n[2]);
		  if (n[0]->num >NB-1)  delete(n[0]->nearnode,n[1]);
		  if (n[1]->num >NB-1) delete(n[1]->nearnode,n[0]);
		  /*
		  for(k=2;k<4;k++)
		    {
		      if (n[k]->num >NB-1)
			{
			  stck=n[k]->nearnode;
			  nx=ny=nsf=0;
			  num=0;
			  while(stck->next!=NULL)
			    {
			      stck=stck->next;
			      nx+=stck->node->x;
			      ny+=stck->node->y;
			      nsf+=stck->node->nsf;
			      num++;
			    }
			  n[k]->x=nx/num;
			  n[k]->y=ny/num;
			  n[k]->nsf=nsf/num;
			}
		    }
		    */	  
		  triangle->node[0]=n[0];
		  triangle->node[1]=n[2];
		  triangle->node[2]=n[3];
		  
		  triangle->neighbour[0]=n20;
		  triangle->neighbour[1]=index;
		  triangle->neighbour[2]=n30;
		  
		  found=0;
		  if (n30!=NULL)
		    {
		      for(k=0;k<3 && !found;k++)
			if (n30->neighbour[k]==index)
			  {
			    found=1;
			    n30->neighbour[k]=triangle;
			  }
		    }
		  
		  index->node[0]=n[1];
		  index->node[1]=n[2];
		  index->node[2]=n[3];
		  
		  index->neighbour[0]=n12;
		  index->neighbour[1]=triangle;
		  index->neighbour[2]=n13;
		  
		  found=0;
		  if (n12!=NULL)
		    {
		      for(k=0;k<3 && !found;k++)
			if (n12->neighbour[k]==triangle)
			  {
			    found=1;
			    n12->neighbour[k]=index;
			  }
		    }
		}
	    }
	}
    }
  return complete;
}


int check_circle (NODE * n[4])
{
  int i;
  float cx,cy,r;
  int flag=0;

  flag=getcircumcircle(n[0],n[1],n[2],&cx,&cy,&r);

  if (!flag)
    if ((n[3]->x-cx)*(n[3]->x-cx)+(n[3]->y-cy)*(n[3]->y-cy) < r) flag=1;

  if (!flag)
   {
     flag=getcircumcircle(n[0],n[1],n[3],&cx,&cy,&r);
     if (!flag)
       if ((n[2]->x-cx)*(n[2]->x-cx)+(n[2]->y-cy)*(n[2]->y-cy) < r) flag=1;
   }

 return flag;
}
/**************************************************************************/
/******* function which gets the circum centre and circum radius   ********/

int getcircumcircle(NODE *N1,NODE *N2,NODE *N3,float *cx,float *cy,float *r)
{
 float x1,y1,x2,y2,x3,y3;
 float xc,yc,r1;
 int flag;

 x1=N1->x;
 y1=N1->y;

 x2=N2->x-N1->x;
 y2=N2->y-N1->y;

 x3=N3->x-N1->x;
 y3=N3->y-N1->y;

 if (x3*y2-x2*y3==0) return 1;
 if (y3*x2-y2*x3==0) return 1;

 xc=0.5*((x3*x3+y3*y3)*y2-(x2*x2+y2*y2)*y3)/(x3*y2-x2*y3);
 yc=0.5*((x3*x3+y3*y3)*x2-(x2*x2+y2*y2)*x3)/(y3*x2-y2*x3);

 r1=xc*xc+yc*yc;
 xc+=x1;
 yc+=y1;
 *cx=xc;
 *cy=yc;
 *r=r1;
 return 0;
}

void swap(TRIANGLE **a,TRIANGLE **b)
{
 TRIANGLE *temp;
 temp=*a;
 *a=*b;
 *b=temp;
}

float distance(NODE *N1,NODE *N2)
{
 return sqrt((N1->x-N2->x)*(N1->x-N2->x)+(N1->y-N2->y)*(N1->y-N2->y));
}


int make_newtriangles1(int num_triangles,TRIANGLE *triangle_head,
		      NODE *node_head)
{
  int i,j,k,k1,found,kk;
  float cx,cy,r,cnsf;
  NODE *node;
  TRIANGLE *triangle,*temp,*maxtri,*temp_tri[3];
  float maxr;
  DELSTACK *stack_head,*stack;

  stack_head=(DELSTACK *)malloc(sizeof(DELSTACK));
  stack_head->next=NULL;
  
  k=num_triangles;
  
  node=node_head;
  while(node->next!=NULL)
    node=node->next;

  triangle=triangle_head;
  
  while(1)
    {
      maxr=0;
      maxtri=NULL;
      triangle=triangle_head;
      while(triangle->next!=NULL)
	{
	  triangle=triangle->next;
	  if (check_nsf(triangle))
	    {
	      getcircumcircle(triangle->node[0],triangle->node[1],
			      triangle->node[2],&cx,&cy,&r);
	      if (r > maxr) 
		{
		  maxr=r;
		  maxtri=triangle;
		}
	    }
	  
	}
      if (maxtri==NULL) break;
      triangle=maxtri;
    

      cx=(triangle->node[0]->x+
	  triangle->node[1]->x+
	  triangle->node[2]->x)/3;
	    
      cy=(triangle->node[0]->y+
	  triangle->node[1]->y+
	  triangle->node[2]->y)/3;
      
      cnsf=0;
      for(kk=0;kk<3;kk++)
	cnsf+=triangle->node[kk]->nsf;
      
      cnsf=cnsf/3.0;
      
      node->next=(NODE *)malloc(sizeof(NODE));
      node=node->next;
      node->next=NULL;
      
      node->x=cx;
      node->y=cy;
      node->nsf=cnsf;
      node->num=number;
      node->nearnode=(STACK *)malloc(sizeof(STACK));
      node->nearnode->next=NULL;
	
      for(k1=0;k1<3;k1++)
	push(node->nearnode,triangle->node[k1]);
      
      for(k1=0;k1<3;k1++)
	if (triangle->node[k1]->num >NB-1)
	  push(triangle->node[k1]->nearnode,node);
      
      number++;
      
      for(i=0;i<3;i++)
	temp_tri[i]=triangle->neighbour[i];
      
      insert_node(triangle,triangle_head,node);
            
      if (temp_tri[0]!=NULL)
	push_tr(stack_head,triangle_head->next,temp_tri[0]);
      if (temp_tri[1]!=NULL)
	push_tr(stack_head,triangle_head->next->next,temp_tri[1]);
      if (temp_tri[2]!=NULL)
	push_tr(stack_head,triangle,temp_tri[2]);
      make_delaunay(stack_head);
      k=k+2;
      //if (k%1000==0) printf("number of triangles=%d\n",k);
    }

  stack=stack_head;
  while(stack_head!=NULL)
    {
      stack_head=stack->next;
      free(stack);
      stack=stack_head;
    }

  num_triangles=k;
  return num_triangles;
}

int check_nsf(TRIANGLE *triangle)
{
 int n[4];
 int k,j;
 TRIANGLE *temp;
 float x1,y1,x2,y2,x3,y3,px,py,cx,cy,nsf1,nsf2,nsf3,ns;

 x1=triangle->node[0]->x;
 x2=triangle->node[1]->x;
 x3=triangle->node[2]->x;
    	    	    
 y1=triangle->node[0]->y;
 y2=triangle->node[1]->y;
 y3=triangle->node[2]->y;

 nsf1=triangle->node[0]->nsf;
 nsf2=triangle->node[1]->nsf;
 nsf3=triangle->node[2]->nsf;

 cx=(x1+x2+x3)*0.33333;
 cy=(y1+y2+y3)*0.33333;
 ns=(nsf1+nsf2+nsf3)/3;
  
 if ((cx-x1)*(cx-x1)+(cy-y1)*(cy-y1)<ns*ns || 
     (cx-x2)*(cx-x2)+(cy-y2)*(cy-y2)<ns*ns || 
     (cx-x3)*(cx-x3)+(cy-y3)*(cy-y3)<ns*ns)
   return 0;
 
 for(j=0;j<3;j++)
   {
     temp=triangle->neighbour[j];
     if (temp!=NULL)
       {
	 px=(temp->node[0]->x
	     +temp->node[1]->x
	     +temp->node[2]->x)/3;
	 py=(temp->node[0]->y+
	     temp->node[1]->y+
	     temp->node[2]->y)/3;
	 if ((px-x1)*(px-x1)+(py-y1)*(py-y1) < ns*ns) return 0;
       }
   }
 return 1;
}

void push(STACK *stack_head,NODE *N)
{
  STACK *temp;
  
  temp=stack_head->next;
  stack_head->next=(STACK *)malloc(sizeof(STACK));
  stack_head->next->node=N;
  stack_head->next->next=temp;
}

void delete(STACK *stack_head,NODE *N)
{
  STACK *temp,*temp1;
  
  temp=stack_head;
  
  while(temp->next!=NULL)
    {
      temp1=temp;
      temp=temp->next;
      if (temp->node==N)
	{
	  temp1->next=temp->next;
	  free (temp);
	  return;
	}
    }
  return;
}


int smooth2(NODE *node_head)
{
 int i,j,count=0,num;
 float sumx,sumy,sumz,sumnsf;
 float EPS1=1e-3,L2;
 float n1,n2;
 STACK *stck;
 NODE *node,*tnode;

 i=0;
 node=node_head;
 while(i<NB)
   {
     node=node->next;
     i++;
   }
 tnode=node;
 do
   {
     L2=0;
     node=tnode;
     while(node->next!=NULL)
       {
	 node=node->next;
	 sumx=sumy=sumz=sumnsf=0;
	 stck=node->nearnode;
	 num=0;
	 while(stck->next!=NULL)
	   {
	     stck=stck->next;
	     sumx+=stck->node->x;
	     sumy+=stck->node->y;
	     sumz+=1/stck->node->nsf;
	     num++;
	   }
	 n1=sumx/num;
	 n2=sumy/num;
	 node->nsf=num/sumz;
	 L2+=((node->x-n1)*(node->x-n1)+(node->y-n2)*(node->y-n2));
	 node->x1=n1;
	 node->y1=n2;
       }
     node=tnode;
     while(node->next !=NULL) 
     {
        node=node->next;
        node->x=node->x1;
        node->y=node->y1;
     } 
     L2=sqrt(L2/(2*(number-NB)));
     count++;
   }
 while(L2>EPS1 && count <20);
 
 if (L2<EPS1) return 0;
 else return 1;
}

int smooth3(NODE *node_head)
{
 int i,j,count=0,num;
 float sumx,sumy,sumz,sumnsf;
 float EPS1=1e-3,L2;
 float n1,n2;
 STACK *stck;
 NODE *node,*tnode;

 i=0;
 node=node_head;
 while(i<NB)
   {
     node=node->next;
     i++;
   }
 tnode=node;
 do
   {
     L2=0;
     node=tnode;
     while(node->next!=NULL)
       {
	 node=node->next;
	 sumx=sumy=sumz=sumnsf=0;
	 stck=node->nearnode;
	 num=0;
	 while(stck->next!=NULL)
	   {
	     stck=stck->next;
	     sumx+=stck->node->x;
	     sumy+=stck->node->y;
	     sumz+=stck->node->nsf;
	     num++;
	   }
	 n1=sumx/num;
	 n2=sumy/num;
	 node->nsf=sumz/num;
	 L2+=((node->x-n1)*(node->x-n1)+(node->y-n2)*(node->y-n2));
	 node->x=n1;
	 node->y=n2;
       }
     L2=sqrt(L2/(2*(number-NB)));
     count++;
   }
 while(L2>EPS1 && count <10);
 
 if (L2<EPS1) return 0;
 else return 1;
}


void print_fast(int n, int number,TRIANGLE *triangle_head,NODE *node_head)
{
 FILE *fp;
 int i,j;
 TRIANGLE *triangle;
 NODE *node;

 fp = fopen("fast.dat","w");

 fprintf(fp,"%d %d 0\n",number,n);
 
 node=node_head->next;
 for(i=0;i<number;i++,node=node->next)
  fprintf(fp,"%f\n",node->x);
 node=node_head->next;
 for(i=0;i<number;i++,node=node->next)
  fprintf(fp,"%f\n",node->y);
 for(i=0;i<number;i++)
  fprintf(fp,"%f\n",1.00);
  
 triangle=triangle_head->next;
 for(i=0;i<n;i++,triangle=triangle->next)
   fprintf(fp,"%d %d %d\n",triangle->node[0]->num+1,triangle->node[1]->num+1,
	  triangle->node[2]->num+1);
 for(i=0;i<n;i++)
   fprintf(fp,"2\n");
 fclose(fp);
}

void insert_node(TRIANGLE *triangle,TRIANGLE *triangle_head,NODE *node)
{
  TRIANGLE *temp;
  int found, k,k1;

  temp=triangle_head->next;
  triangle_head->next=(TRIANGLE *)malloc(sizeof(TRIANGLE));
  triangle_head->next->next=(TRIANGLE *)malloc(sizeof(TRIANGLE));
  triangle_head->next->node[0]=triangle->node[0];
  triangle_head->next->node[1]=triangle->node[1];
  triangle_head->next->node[2]=node;
  triangle_head->next->neighbour[0]=triangle->neighbour[0];
  triangle_head->next->neighbour[1]=triangle_head->next->next;
  triangle_head->next->neighbour[2]=triangle;
        
  found=0;
  if (triangle->neighbour[0]!=NULL)
    {
      for(k1=0;k1<3 && !found;k1++)
	if (triangle->neighbour[0]->neighbour[k1]==triangle)
	  {
	    found=1;
	    triangle->neighbour[0]->neighbour[k1]=triangle_head->next;
	  }
    }
  
  triangle_head->next->next->node[0]=triangle->node[1];
  triangle_head->next->next->node[1]=triangle->node[2];
  triangle_head->next->next->node[2]=node;
  
  triangle_head->next->next->neighbour[0]=triangle->neighbour[1];
  triangle_head->next->next->neighbour[1]=triangle;
  triangle_head->next->next->neighbour[2]=triangle_head->next;
  
  found=0;
  if (triangle->neighbour[1]!=NULL)
    {
      for(k1=0;k1<3 && !found;k1++)
	if (triangle->neighbour[1]->neighbour[k1]==triangle)
	  {
	    found=1;
	    triangle->neighbour[1]->neighbour[k1]
                    =triangle_head->next->next;
	  }
    }
  
  triangle->node[1]=node;
  triangle->neighbour[0]=triangle_head->next;
  triangle->neighbour[1]=triangle_head->next->next;
  triangle_head->next->next->next=temp;
}  


/** Doing Lawsons edge swap to get the C-surface correct */

void make_delaunay(DELSTACK *stack_head)
{
  TRIANGLE *tri,*tri1; 
  TRIANGLE *n13,*n30,*n20,*n12;
  NODE *n[4];
  int i,k,found;

  while(stack_head->next!=NULL)
    {
     pop_tr(stack_head,&tri,&tri1);
     
     for(i=0;i<3;i++) 
       if (tri->neighbour[i]==tri1) 
	 {
	   n[2]=tri->node[(i+2)%3];
	   n[0]=tri->node[i];
	   n[1]=tri->node[(i+1)%3];
	   
	   n12=tri->neighbour[(i+1)%3];
	   n20=tri->neighbour[(i+2)%3];
	 }

     for(i=0;i<3;i++)
       if (tri1->neighbour[i]==tri)
	 {
	   n[3]=tri1->node[(i+2)%3];
	   
	   n13=tri1->neighbour[(i+1)%3];
	   n30=tri1->neighbour[(i+2)%3];

	   if (tri1->node[i]==n[1])
	     swap(&n13,&n30);
	 }

     if (check_circle(n))
       {  
	 tri->node[0]=n[0];
	 tri->node[1]=n[3];
	 tri->node[2]=n[2];
		  
         tri->neighbour[0]=n30;
	 tri->neighbour[1]=tri1;
	 tri->neighbour[2]=n20;
	 
	 if (n[2]->num >NB-1)  push(n[2]->nearnode,n[3]);
	 if (n[3]->num >NB-1)  push(n[3]->nearnode,n[2]);
	 if (n[0]->num >NB-1)  delete(n[0]->nearnode,n[1]);
	 if (n[1]->num >NB-1) delete(n[1]->nearnode,n[0]);
	 
	 if (n30!=NULL)
	   push_tr(stack_head,tri,n30);
		  
	 found=0;
	 if (n30!=NULL)
	   {
	     for(k=0;k<3 && !found;k++)
	       if (n30->neighbour[k]==tri1)
		 {
		   found=1;
		   n30->neighbour[k]=tri;
		 }
	   }
		  
	 tri1->node[0]=n[1];
	 tri1->node[1]=n[2];
	 tri1->node[2]=n[3];
	 
	 tri1->neighbour[0]=n12;
	 tri1->neighbour[1]=tri;
	 tri1->neighbour[2]=n13;
	 
	 if (n13!=NULL)
	   push_tr(stack_head,tri1,n13);
	 
	 found=0;
	 if (n12!=NULL)
	   {
	     for(k=0;k<3 && !found;k++)
	       if (n12->neighbour[k]==tri)
		 {
		   found=1;
		   n12->neighbour[k]=tri1;
		 }
	   }
       }
    }
}

/** standard push and pop **/

void push_tr(DELSTACK *stack_head,TRIANGLE *tri,TRIANGLE *tri1)
{
  DELSTACK *temp;

  if (tri1!=NULL)
    {
      temp=stack_head->next;
      stack_head->next=(DELSTACK *)malloc(sizeof(DELSTACK));
      stack_head->next->tri=tri;
      stack_head->next->tri1=tri1;
      stack_head->next->next=temp;
    }
}

void pop_tr(DELSTACK *stack_head,TRIANGLE **tri,TRIANGLE **tri1)
{
  DELSTACK *temp;
  
  temp=stack_head->next;
  stack_head->next=temp->next;
  
  *tri1=temp->tri1;
  *tri=temp->tri;
  free(temp);
}
  








