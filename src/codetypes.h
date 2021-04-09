//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef CODETYPES_H
#define CODETYPES_H

//#define MPICH_SKIP_MPICXX
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include "mpi.h"
/*====================================================================*/
/*  Floating point definition                                         */
/*====================================================================*/
# define REAL double

/*====================================================================*/
/*  Define arithmetic constants                                       */
/*====================================================================*/
#define ZERO               0.0e+00
#define ONE                1.0e+00
#define TWO                2.0e+00
#define THREE              3.0e+00
#define FOUR               4.0e+00
#define HALF               0.5e+00
#define THIRD              0.333333333e+00
#define FIFTH              0.2
#define PI                 3.1415926535897932e+00
#define RAD2DEG            (180.0/PI)
#define DEG2RAD            (PI/180.0)
#define BIGVALUE           1.0e+15
#define BIGINT             2147483647
#define TOL                1.0e-10
#define HOLEMAPSIZE        192
/*==================================================================*/
/* inline debugging tools                                             */
/*==================================================================*/
# define TRACEI(x)  printf("#dgsand:\t"#x" =%d\n",x);
# define TRACED(x)  printf("#dgsand:\t"#x" =%.16e\n",x);
# define DGSAND_MIN(x,y)  (x) < (y) ? (x) : (y)
# define DGSAND_MAX(x,y)  (x) > (y) ? (x) : (y)
# define DGSAND_FREE(a1)  {free(a1);a1=NULL;}
# define debug(x,y)  printf("#dgsand:\t"#x"=%d,"#y"=%d\n",x,y);
# define stdwrite(x) if (myid==0) printf("#dgsand:\t"#x"\n");
# define dstr(x) printf("#dgsand:\t"#x"\n");
//# define ditch(x,y) {dstr(x); TRACEI(y); MPI_Abort(MPI_COMM_WORLD,ierr);}
/*====================================================================*/
/*  Numerical Tools                                                   */
/*====================================================================*/
#define DGSAND_Sign(a1,a2) (((a2) < ZERO)? - fabs(a1): fabs(a1))
#define DGSAND_Max(a1,a2) (((a1) >= (a2))? (a1): (a2))
#define DGSAND_Min(a1,a2) (((a1) <= (a2))? (a1): (a2))
#define DGSAND_Abs(aa) (((aa) >= 0)? (aa): -(aa))
#define DGSAND_Round(aa) (int) ((fabs((aa) - floor(aa)) >= HALF)? ceil(aa): floor(aa))
#define DGSAND_swap(a,b) { a=a+b;b=a-b;a=a-b;}
/*********************************************************************/
/* Code specific types */
/*********************************************************************/
typedef struct ELEMENT
{
 int elementType;  // type of element
 int pdegree;      // p degree
 int nbasis;       // number of basis functions   
 int nquad;        // number of quadrature points (volume)
 int nfluxquad;    // number of flux quadrature points (surface)
 int *cell2node;   // cell 2 node connectivity
 int *cellface ;   // cell 2 face connectivity
 double *a;        // modal coefficients (unknowns) a(nbasis,nfield)
 double *b;        // basis values at quadrature points b(nbasis,nquad)
 double *bf;       // basis values at flux points (nbasis,nfluxquad)
 double *gradb     // gradient of basis at quadrature points gradb(nbasis,nquad)
 double *J;        // Jacobian values at quadrature points J(nquad,d,d)
 double *detJ;     // determinant of the Jacobian detJ(nquad)
 double *M;        // mass matrix (nbasis,nbasis)
 double *residual; // residual (nbasis,nfield)
} ELEMENT;

typedef struct FACE
{
 int *face2node // face connectivity to nodes
 int *elem;     // elements on either side of the face
 int *eflux;    // element flux on either side eflux(2,nfield)
 int *rflux;    // riemann flux (nfield)
} FACE;

typedef struct REFELEMENT
{
 double *b;
 double *bf;
 double *gradb;
} REFLEMENT;


typedef struct GRID
{
 int nnodes,nelem,nfaces; // number of nodes, elements and faces
 double *x;               // coordinate list
 int *cell2node;          // cell2node connectivity   
 int *cell2face;          // cell2face connectivity
 int *face;               // face connectivity 
 double *xb;             // coordinate basis coefficients
 double *b;              // ref element basis value b(nref,nbasis)
 double *bf;             // flux basis value
 double *gradb;          // gradient of basis function
 double *J;              // Jacobians for all elements + all quadrature J(nelem,nquad,d,d)
 double *detJ;           // determinant of Jacobian detJ(nelem,nquad)         
 double *M;              // Mass matrix for each element (nelem,nbasis,nbasis)
} GRID;


typedef struct SOLUTION
{
 double *a;                  // solution coefficients
 double *eflux;              // face flux eflux(nelem,2,nfield)
 double *rflux;              // riemann flux rflux(nelem,nfield)
 double *residual;           // residual (nelem,nfield)
} SOLUTION;

typedef struct FLOWPARAM
{
 double V[3];                // free stream velocity
 double dt;                  // time step
 double CFL;                 // CFL
 int initcond;               // initial condition
} FLOWPARAM;
 
