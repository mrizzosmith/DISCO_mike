#include "../mike.h"

static double rate = 0.0;
static double dexp = 0.0;
static double q = 0.0;
static double t_max = 0.0;
static double mach = 0.0;
static double eps_frac = 0.0;
static double r_init = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2;

   rate  = theDomain->theParList.Drift_Rate;
   dexp  = theDomain->theParList.Drift_Exp;
   q     = theDomain->theParList.Mass_Ratio;
   t_max = theDomain->theParList.t_max;
   mach  = theDomain->theParList.Disk_Mach;
   eps_frac = theDomain->theParList.eps_frac;
   r_init = theDomain->theParList.r_init;

}

int planet_motion_analytic( void ){
   return(1);
}

double drift_pos( double R , double t ){
// MRS variable r_init corresponds to setting final separation of drift (double check before using)
   return( r_init * pow( 1. + R*( t - t_max )/dexp , dexp ) );
}

void initializePlanets( struct planet * thePlanets ){

   double r = drift_pos( rate , 0.0 );

   double mu = q/(1+q);

   thePlanets[0].M     = 1.0 - mu; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = pow(r, -1.5); 
   thePlanets[0].r     = r*mu; 
   thePlanets[0].phi   = M_PI; 
   thePlanets[0].eps   = 0.0;

   thePlanets[1].M     = mu; 
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = pow(r,-1.5); 
   thePlanets[1].r     = r*(1-mu); 
   thePlanets[1].phi   = 0.0; 
   //amd: smoothing length is some fraction of scale height
   thePlanets[1].eps   = eps_frac/mach;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){

   double r = drift_pos( rate , t+dt );

   double mu = q/(1+q);

   thePlanets[0].phi  += thePlanets[0].omega*dt;
   thePlanets[0].r     = r*mu;

   thePlanets[1].r     = r*(1-mu);
   thePlanets[1].omega = pow(r,-1.5);
   thePlanets[1].phi  += thePlanets[1].omega*dt;
   thePlanets[1].eps   = eps_frac * r/mach;

}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

