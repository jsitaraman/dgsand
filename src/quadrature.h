

/*  gauss point locations                                                       */
/*         (r,s,w) in UCS [0,1] coordinates system for triangles                */
/*         TODO these weights need to be halved.. currently fixed in code       */
/*                                r                 s              weight       */

double tri0[3]  = { 3.333333333333335e-01, 3.333333333333335e-01, 0.500000000000000e+00 };

double tri1[12] ={  1.666666666666665e-01, 1.666666666666665e-01, 1.666666666666666e-01,
		    1.666666666666665e-01, 6.666666666666665e-01, 1.666666666666667e-01,
		    6.666666666666665e-01, 1.666666666666665e-01, 1.666666666666667e-01 };

double tri2[18] ={  4.459484909159650e-01, 4.459484909159650e-01, 1.116907948390056e-01,
		    4.459484909159650e-01, 1.081030181680700e-01, 1.116907948390057e-01,
		    1.081030181680700e-01, 4.459484909159650e-01, 1.116907948390057e-01,
		    9.157621350977102e-02, 9.157621350977102e-02, 5.497587182766102e-02,
		    9.157621350977102e-02, 8.168475729804590e-01, 5.497587182766102e-02,
		    8.168475729804590e-01, 9.157621350977102e-02, 5.497587182766103e-02};



/* (r,w) of GL points in 1D - weights are correct */
double gl0[2] = {  5.000000000000000e-01, 1.000000000000000e+00 };

double gl1[4] = {  2.113248654051870e-01, 5.000000000000000e-01,
		   7.886751345948130e-01, 5.000000000000000e-01  };

double gl2[6]= {  1.127016653792585e-01, 2.777777777777780e-01,
		  5.000000000000000e-01, 4.444444444444445e-01,
		  8.872983346207415e-01, 2.777777777777780e-01 };

double *gauss[1][3]={tri0,tri1,tri2};
double *gaussgl[1][3]={gl0,gl1,gl2};

/* number of quadrature points per element type and polynomial order  */
int ngElem[1][5]={1,3,6,6,6};
int ngGL[1][5]={1,2,3,3,3};
/* map of polynomial order to gauss-quadrature type for volume integration */
int p2g[1][5]={0,1,2,2,2};
/* map of polynomial order to gauss-quadrature type for face integration */
int p2gf[1][5]={0,1,2,2,2};
