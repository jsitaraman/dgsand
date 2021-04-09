/*
hierarchic basis from 
https://www.math.u-bordeaux.fr/~durufle/montjoie/triangle.php#TriangleHierarchic
*/

#define jac2(z) (3*(1+2.5*((z)-1)+1.25*((z)-1)*((z)-1)))
#define jac2d(z) (7.5+7.5*((z)-1))
  
/* hierarchical basis for triangles upto p=4 */

double tri_h0(double *u)  {double r=u[0],s=u[1]; return 1-r-s;}
double tri_h1(double *u)  {double r=u[0],s=u[1]; return     r;}
double tri_h2(double *u)  {double r=u[0],s=u[1]; return     s;}
double tri_e0(double *u)  {double r=u[0],s=u[1]; return r*(1-r-s);}
double tri_e1(double *u)  {double r=u[0],s=u[1]; return r*s;}
double tri_e2(double *u)  {double r=u[0],s=u[1]; return s*(1-r-s);}
double tri_e3(double *u)  {double r=u[0],s=u[1]; return (1-r-s)*r*(s-1)*2;}
double tri_e4(double *u)  {double r=u[0],s=u[1]; return r*s*(s-r)*2;}
double tri_e5(double *u)  {double r=u[0],s=u[1]; return (1-r-s)*s*(r-1)*2;}
double tri_e6(double *u)  {double r=u[0],s=u[1]; return (1-r-s)*r*jac2(s-1);}
double tri_e7(double *u)  {double r=u[0],s=u[1]; return r*s*jac2(s-r);}
double tri_e8(double *u)  {double r=u[0],s=u[1]; return s*(1-r-s)*jac2(r-1);}
double tri_c00(double *u) {double r=u[0],s=u[1]; return r*s*(1-r-s);}
double tri_c10(double *u) {double r=u[0],s=u[1]; return r*s*(1-r-s)*2*(2*r-1+s);}
double tri_c01(double *u) {double r=u[0],s=u[1]; return r*s*(1-r-s)*2*(2*s-1);}

/* derivatives of basis functions w.r.t to "r" and "s" */

double tri_h0r(double *u)  {return -1;}
double tri_h0s(double *u)  {return -1;}

double tri_h1r(double *u)  {return 1; }
double tri_h1s(double *u)  {return 0; }

double tri_h2r(double *u) { return 0;}
double tri_h2s(double *u) { return 1;}

double tri_e0r(double *u)  {double r=u[0],s=u[1]; return (1.0-2*r-s); }
double tri_e0s(double *u)  {double r=u[0],s=u[1]; return  -r;}

double tri_e1r(double *u)  {double r=u[0],s=u[1]; return s;}
double tri_e1s(double *u)  {double r=u[0],s=u[1]; return r;}

double tri_e2r(double *u)  {double r=u[0],s=u[1]; return -s;}
double tri_e2s(double *u)  {double r=u[0],s=u[1]; return (1.0-r-2*s);}

double tri_e3r(double *u)  {double r=u[0],s=u[1]; return 2*(s-1)*(1-2*r-s);}
double tri_e3s(double *u)  {double r=u[0],s=u[1]; return 2*2*r*(2-r-2*s);}

double tri_e4r(double *u)  {double r=u[0],s=u[1]; return 2*s*(s-2*r);}
double tri_e4s(double *u)  {double r=u[0],s=u[1]; return 2*(r*(2*s-r));}

double tri_e5r(double *u)  {double r=u[0],s=u[1]; return 2*(s*(2-2*r-s));}
double tri_e5s(double *u)  {double r=u[0],s=u[1]; return 2*((r-1)*(1-r-2*s));}

double tri_e6r(double *u)  {double r=u[0],s=u[1]; return (1-2*r-s)*jac2(s-1);}
double tri_e6s(double *u)  {double r=u[0],s=u[1]; return r*(-jac2(s-1)+(1-r-s)*jac2d(s-1));}

double tri_e7r(double *u)  {double r=u[0],s=u[1]; return s*(jac2(s-r)-r*jac2d(r-1));}
double tri_e7s(double *u)  {double r=u[0],s=u[1]; return r*(jac2(s-r)+s*jac2d(s-r));}

double tri_e8r(double *u)  {double r=u[0],s=u[1]; return s*(-jac2(r-1)+(1-r-s)*jac2d(r-1));}
double tri_e8s(double *u)  {double r=u[0],s=u[1]; return (1-r-2*s)*jac2(r-1);}

double tri_c00r(double *u)  {double r=u[0],s=u[1]; return s*(1-2*r-s);}
double tri_c00s(double *u)  {double r=u[0],s=u[1]; return r*(1-r-2*s);}

double tri_c10r(double *u)  {double r=u[0],s=u[1]; return 2*s*((1-r-s)*(2*r-1+s)-r*(2*r-1+s)+2*r*(1-r-s));}
double tri_c10s(double *u)  {double r=u[0],s=u[1]; return 2*r*((1-r-s)*(2*r-1+s)-s*(2*r-1+s)+2*s*(1-r-s));}

double tri_c01r(double *u)  {double r=u[0],s=u[1]; return 2*s*(2*s-1)*(1-2*r-s);}
double tri_c01s(double *u)  {double r=u[0],s=u[1]; return 2*r*((1-r-s)*(2*s-1)-s*(2*s-1)+s*(1-r-s)*2);}

/* basis function array */
base tribasis[15]={&tri_h0,&tri_h1,&tri_h2,
		   &tri_e0,&tri_e1,&tri_e2,
		   &tri_e3,&tri_e4,&tri_e5,
		   &tri_c00,
		   &tri_e6,&tri_e7,&tri_e8,
		   &tri_c01,&tri_c10};

/* basis derivative function array */
base tribasis_d[30]={&tri_h0r,&tri_h0s,
		     &tri_h1r,&tri_h1s,
		     &tri_h2r,&tri_h2s,
		     &tri_e0r,&tri_e0s,
		     &tri_e1r,&tri_e1s,
		     &tri_e2r,&tri_e2s,
		     &tri_e3r,&tri_e3s,
		     &tri_e4r,&tri_e4s,
		     &tri_e5r,&tri_e5s,
		     &tri_c00r,&tri_c00s,
		     &tri_e6r,&tri_e6s,
		     &tri_e7r,&tri_e7s,
		     &tri_e8r,&tri_e8s,
		     &tri_c10r,&tri_c10s,
		     &tri_c01r,&tri_c01s};

/* number of basis for each order up to p=4 */
int order2basis_tri[5]={1,3,6,10,15};


/* location of nodes in triangle grids */
/* r,s values of physical coordinate locations for elements             */
/* p = 0 */
double tri0loc[2]={ 0.333333333333333, 0.33333333333333};

/* p = 1 */
double tri1loc[6] = { 0.0, 0.0,
               1.0, 0.0,
	       0.0, 1.0} ;

/* p = 2 */
double tri2loc[12] = { 0.0, 0.0,
		1.0, 0.0,
		0.0, 1.0,
		0.5, 0.0,
		0.5, 0.5,
	        0.0, 0.5};

/* p = 3 */
double tri3loc[20] = { 0.0, 0.0,
		1.0, 0.0,
		0.0, 1.0,
		0.33333333333333333, 0.0,
		0.66666666666666667, 0.3333333333333333,	
		0.0, 0.66666666666666667,
		0.66666666666666667, 0.0,
	        0.33333333333333333, 0.6666666666666667,
		0.0, 0.33333333333333333,
		0.33333333333333333, 0.33333333333333333,
                };
/* p = 4 */
double tri4loc[30] = { 0.0, 0.0,
		1.0, 0.0,
		0.0, 1.0,
		0.25, 0.0,
		0.75, 0.25,	
		0.0,  0.75,
		0.5, 0.0,
	        0.5, 0.5,
		0.0,  0.5,
		0.25, 0.25,
		0.75, 0.0,
                0.25, 0.75,
                0.0,  0.25,
		0.5,  0.25,
		0.25, 0.5};
                
/* element coord locations */
double *eloctri[5]={tri0loc,tri1loc,tri2loc,tri3loc,tri4loc};

/* tranforms between face coordinate system and element coordinate system */
/* these are dr/de ds/de where e is the edge coordinate direction for each face*/
 
double face2tri[6]={1,0,-1,1,0,-1};
//double face2tri[6]={0,-1,1,1,-1,0};
