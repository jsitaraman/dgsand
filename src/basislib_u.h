/* simple definition of a basis function, takes in an array of barycentric coordinates as argument */
/* should be exactly same for any order and any dimension */

typedef double (* base)(double *);

/* hierarchical triangle basis, TODO, add more elements here */
#include "tribasis.h"

/* combine basis and basis derivative functions
   from all element libraries   */

base *basis[1]={tribasis};
base *basis_d[1]={tribasis_d};

/* location of control nodes for each element type grid */
double **eloc[1]={eloctri};

/* mapping between order and number of basis for each element type */
int *order2basis[1]={order2basis_tri};

/* faces per element */
int facePerElem[1]={3};

/* face transforms per element */
double *face2elem[1]={face2tri};
