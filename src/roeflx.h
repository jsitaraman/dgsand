//***************************************************************************
//
// Roe's approximate Riemann solver
//
// FLUX = 0.5*((FL+FR)-|A|*(QR-QL))
//
// This implementation is from
// 
// Vatsa, V. N., Thomas, J., L. and Wedan, B.,W.,"Navier-Stokes Computations 
// of Prolate Spheroids at Angle of Attack," AIAA 87-2627
//
//
// Note:
//
// Includes grid speed terms
// Closely follows implementation in UMTURNS code
// The left and right states are conservative variables
//
// T|lambda|Tinv is abbreviated for decreasing floating point op count
//
//J. Sitaraman 06/09/09 (F90 version)
// J. Sitaraman 02/27/11 (C version)
//
// ***************************************************************************
#define roe_max(a,b) (a>b)?a:b
void roeflx(double *leftState,double *rightState, double *flux,double *ds,double faceVel)
{  
  int irho=0,irhou=1,irhov=2,ie=3;
  double  eps,rlft,ulft,vlft,plft;
  double  rlfti,rulft,rvlft,uvl,elft,hlft,clft;
  double  rrht,urht,vrht,prht;
  double  rrhti,rurht,rvrht,uvr,erht,hrht,crht;
  double  rat,rati,rav,uav,vav,hav,uv,cav;
  double  aq1,aq2,aq3,aq4,ri1,ri2,ri3,rr2,rr,r0,r1,r2,r3;
  double  uu,c2,c2i,auu,aupc,aumc,uulft,uurht,upclft,upcrht;
  double  umclft,umcrht,dauu,dauus,daupc,daumc,daumcs,rcav,aquu;
  double  daupcs,c2ih,ruuav,b1,b2,b3,b4,b5,b6,b7,aj;
  double  plar,eplft,eprht,fssub;
  double  uLi[4],uRi[4],nx,ny,area;
  double gm1=gamma-1;
  double specRadius;
  //
  eps=1e-14;
  //  
  rlft = leftState[irho];
  ulft = leftState[irhou]/leftState[irho];
  vlft = leftState[irhov]/leftState[irho];
  plft = gm1*(leftState[ie]-0.5*rlft*(ulft*ulft+vlft*vlft));
  //
  rlfti = 1.0/rlft;
  rulft = rlft*ulft;
  rvlft = rlft*vlft;
  uvl = 0.5*( ulft*ulft + vlft*vlft );
  elft = plft/gm1 + rlft*uvl;
  hlft = ( elft + plft )*rlfti;
  clft = sqrt( gm1*( hlft - uvl ) );
  //
  rrht = rightState[irho];
  urht = rightState[irhou]/rightState[irho];
  vrht = rightState[irhov]/rightState[irho];
  prht = gm1*(rightState[ie]-0.5*rrht*(urht*urht+vrht*vrht));
//  printf("Inside Roe flux\n L state %f %f %f %f\n R state %f %f %f %f\n",rlft,ulft,vlft,plft,rrht,urht,vrht,prht); 
  rrhti = 1.0/rrht;
  rurht = rrht*urht;
  rvrht = rrht*vrht;
  uvr = 0.5*( urht*urht + vrht*vrht );
  erht = prht/gm1 + rrht*uvr;
  hrht = ( erht + prht )*rrhti;
  crht = sqrt( gm1*( hrht - uvr ) );
  //
  rat  = sqrt( rrht*rlfti );
  rati = 1.0/( rat + 1. );
  rav  =   rat*rlft;
  uav  = ( rat*urht + ulft )*rati;
  vav  = ( rat*vrht + vlft )*rati;
  hav  = ( rat*hrht + hlft )*rati;
  uv   = 0.5*( uav*uav + vav*vav );
  cav  = sqrt( gm1*( hav - uv ) );
  //
  aq1  = rrht - rlft;
  aq2  = urht - ulft;
  aq3  = vrht - vlft;
  aq4  = prht - plft;
  //
  ri1 = ds[0];
  ri2 = ds[1];
  ri3 = faceVel;
  rr2 = ri1*ri1 + ri2*ri2;
  rr  = sqrt( rr2 );
  r0  = 1.0 / rr;
  r1  = ri1*r0;
  r2  = ri2*r0;
  r3  = ri3*r0;
  //
  uu  = r1*uav + r2*vav + r3;
  c2  = cav*cav;
  c2i = 1.0/c2;
  //
  auu   = fabs( uu    );
  aupc  = fabs( uu+cav );
  aumc  = fabs( uu-cav );
  //
  uulft = r1*ulft + r2*vlft + r3;
  uurht = r1*urht + r2*vrht + r3;
  upclft= uulft + clft;
  upcrht= uurht + crht;
  umclft= uulft - clft;
  umcrht= uurht - crht;
  //
  // entropy fix to minimize expansion shocks
  //
  dauu = 4.*(uurht-uulft)+eps;
  dauus = roe_max(dauu,0.0);
  if( auu <= 0.5*dauus ) auu = auu*auu/dauu+0.25*dauu;
  //
  daupc = 4.*(upcrht-upclft)+eps;
  daupcs = roe_max(daupc,0.0);
  if( aupc <= 0.5*daupcs ) aupc = aupc*aupc/daupc+0.25*daupc;
  //
  daumc = 4.*(umcrht-umclft)+eps;
  daumcs = roe_max(daumc,0.0);
  if( aumc <= 5*daumcs ) aumc = aumc*aumc/daumc+0.25*daumc;
  //
  specRadius=(auu+cav)*rr;
  //
  rcav = rav*cav;
  aquu = uurht - uulft;
  c2ih = 0.5*c2i;
  ruuav= auu*rav;
  b1   = auu*( aq1 - c2i*aq4 );
  b2   = c2ih*aupc*( aq4 + rcav*aquu );
  b3   = c2ih*aumc*( aq4 - rcav*aquu );
  b4   = b1 + b2 + b3;
  b5   = cav*( b2 - b3 );
  b6   = ruuav*( aq2 - r1*aquu );
  b7   = ruuav*( aq3 - r2*aquu );
  //
  aq1 = b4;
  aq2 = uav*b4 + r1*b5 + b6;
  aq3 = vav*b4 + r2*b5 + b7;
  aq4 = hav*b4 + ( uu-r3 )*b5 + uav*b6 + vav*b7 - c2*b1/gm1;
  //
  aj    = 0.5*rr;
  plar  = plft + prht;
  eplft = elft + plft;
  eprht = erht + prht;
  //
  // flux = 0.5*((fL+fR) + |A|(uR-uL)).ds
  //
  //
//  printf("aq1 = %f %f %f %f\n", aq1,aq2,aq3,aq4);
//  aq1=aq2=aq3=aq4=0;
  flux[0] = aj*( rlft*uulft+rrht*uurht-aq1 );
  flux[1] = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 );
  flux[2] = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 );
  flux[3] = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 );

//printf("fluxes = %f %f %f %f\n",flux[0],flux[1],flux[2],flux[3]); 
}


