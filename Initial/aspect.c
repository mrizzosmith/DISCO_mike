
#include "../mike.h"

static double gam  = 0.0;
static double alpha   = 0.0;
static double Mach = 0.0;
static double rho_floor = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   alpha   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   rho_floor = theDomain->theParList.Density_Floor;

}

double get_cs2( double ); 

void initial( double * prim , double * x ){

   double r = x[0];

   double cs2 = get_cs2( r );

   double rho = pow(r,-0.5);

   if( rho < rho_floor ) rho = rho_floor;

   double Pp = rho*cs2;
   double omega02 = 1.0/pow(r,3.);
   double omegaP2 = 1.5*cs2/r/r;

   double omega = sqrt( omega02 - omegaP2 );

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*alpha*cs2/omega/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
