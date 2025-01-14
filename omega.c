
#include "mike.h"

static double Mach = 0.0;

void setDiskParams( struct domain * theDomain ){
   Mach = theDomain->theParList.Disk_Mach;
}

double mesh_om( double r ){
//   double n = 8.0;
//   double omega = 1./pow( pow( r , 1.5*n ) + 1. , 1./n );
   double omega = pow(r,-1.5);
//   double omega = 1.0;
   return( omega );
}

double get_om( double r ){
   return(0.0);
   return( 1.0/pow(r,1.5) );
//   return( 10.*exp(-.5*r*r) );
}
  
double get_om1( double r ){
   return(0.0);
   return( -1.5/pow(r,2.5) );
//   return( -10.*r*exp(-.5*r*r) );
}

double get_cs2( double r ){
//   double nu = .5;
//   MRS This cs2 corresponds to aspect disk
   return( 1./Mach/Mach/r );

//   MRS This cs2 corresponds to a constant nu viscosity (but double check / derive this)
//   return( 1./Mach/Mach );
	
//   return(1.0);
}

