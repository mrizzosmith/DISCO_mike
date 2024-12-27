
#include "../mike.h"

static double q_planet = 1.0;
static double Mach = 1.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2; 
   q_planet = theDomain->theParList.Mass_Ratio;
   Mach = theDomain->theParList.Disk_Mach;

}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){
   
   // r is total separation
   double r = 1.0;
   double mu = q_planet/(1.0 + q_planet);
   double om = pow(r,-1.5);
   
   thePlanets[0].M     = (1.0 - mu);
   thePlanets[0].vr    = 0.0;
   thePlanets[0].omega = om;
   thePlanets[0].r     = r*mu;
   thePlanets[0].phi   = M_PI;
   thePlanets[0].eps   = 0.0;

   thePlanets[1].M     = 0.0;
   thePlanets[1].vr    = 0.0;
   thePlanets[1].omega = om;
   thePlanets[1].r     = r*(1.0 - mu);
   thePlanets[1].phi   = 0.0;
   thePlanets[1].eps   = 0.5/Mach;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[1].phi += thePlanets[1].omega*dt;
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

