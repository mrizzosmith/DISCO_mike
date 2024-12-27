
#include "../mike.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   double r   = x[0];
   double phi = x[1];

   double rho  = 1.0;
   double Pp   = 0.1;
   double phi0 = M_PI/4.;

   double b = 1.0;

   if( r<0.1 ){
      Pp  = 10.0;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0;

   prim[BRR] =  b*cos(phi-phi0);
   prim[BPP] = -b*sin(phi-phi0);
   prim[BZZ] = 0.0;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
