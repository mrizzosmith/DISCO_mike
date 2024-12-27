#include "mike.h"

static double t_sink_factor;
static double mach;
static double alpha;

void setGenericParams( struct domain * theDomain ){
  t_sink_factor  = theDomain->theParList.t_sink_factor;
  mach  = theDomain->theParList.Disk_Mach;
  alpha = theDomain->theParList.viscosity;
}

double PHI_ORDER = 2.0;

// amd
double get_dp( double , double );

double phigrav( double M , double r , double eps ){
   double n = PHI_ORDER;
   return( M/pow( pow(r,n) + pow(eps,n) , 1./n ) ) ;
   // amd: PW potential!
   // scaling for rS depends on post-processing scale factor
   // here r0 = 8rS, so 1rS = 0.125 r0 for M=1
   //return( M/(pow( pow(r,n) + pow(eps,n) , 1./n ) - 0.125*M) ) ;
   // Here as a test r0 = 16rS, so 1rS = 0.0625 for M=1
   //return( M/(pow( pow(r,n) + pow(eps,n) , 1./n ) - 0.0625*M) ) ;
}

double fgrav( double M , double r , double eps ){
   double n = PHI_ORDER;
   return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
   // amd: Changing to Paczynski-Wiita grav force! d/dr of potential
   // Taking out that phi order option
   // pow pow pow pow pow 
   // here r0 = 8rS, so 1rS = 0.125 r0 for M=1
   //return( M*r/(pow(pow(r,2.) + pow(eps,2.) ,0.5) *pow(pow(pow(r,2.) + pow(eps,2.) ,0.5) - 0.125*M, 2.)));
   // Here as a test r0 = 16rS, so 1rS = 0.0625 for M=1
   //return( M*r/(pow(pow(r,2.) + pow(eps,2.) ,0.5) *pow(pow(pow(r,2.) + pow(eps,2.) ,0.5) - 0.0625*M, 2.)));
}

void adjust_gas( struct planet * pl , double * x , double * prim , double gam ){

   double r   = x[0];
   double phi = x[1];

   double rp = pl->r;
   double pp = pl->phi;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy);

   double pot = phigrav( pl->M , script_r , pl->eps );

   double c2 = gam*prim[PPP]/prim[RHO];
   double factor = 1. + (gam-1.)*pot/c2;

   prim[RHO] *= factor;
   prim[PPP] *= pow( factor , gam );

}

void planetaryForce( struct planet * pl , double r , double phi , double z , double * fr , double * fp , double * fz , int mode ){

   z = 0.0;

   double rp = pl->r;
   double pp = pl->phi;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy+z*z);
   double script_r_perp = sqrt(dx*dx+dy*dy);

   double f1 = -fgrav( pl->M , script_r , pl->eps );

   // geoff found a bug!
   double cosa = dx/script_r_perp;
   double sina = dy/script_r_perp;

   double cosap = cosa*cosp+sina*sinp;
   double sinap = sina*cosp-cosa*sinp;

   if( mode==1 ){
      cosap = cosa*cos(pp)+sina*sin(pp);
      sinap = sina*cos(pp)-cosa*sin(pp);
   }
/*
   double rH = rp*pow( pl->M/3.,1./3.);
   double pd = 0.8; 
   double fd = 1./(1.+exp(-( script_r/rH-pd)/(pd/10.)));
*/

   double sint = script_r_perp/script_r;
   double cost = z/script_r;

   *fr = cosap*f1*sint; //*fd;
   *fp = sinap*f1*sint; //*fd;
   *fz = f1*cost;

}

void planet_src( struct planet * pl , double * prim , double * cons , double * xp , double * xm , double dVdt ){

   double rp = xp[0];
   double rm = xm[0];
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double omega = prim[UPP];
   
   double r = 0.5*(rp+rm);
   double vp  = r*omega;
   double dphi = get_dp(xp[1],xm[1]);
   double phi = xm[1] + 0.5*dphi;
   double z = .5*(xp[2]+xm[2]);

   double Fr,Fp,Fz;
   planetaryForce( pl , r , phi , z , &Fr , &Fp , &Fz , 0 );

   cons[SRR] += rho*Fr*dVdt;
   cons[SZZ] += rho*Fz*dVdt;
   cons[LLL] += rho*Fp*r*dVdt;
   cons[TAU] += rho*( Fr*vr + Fz*vz + Fp*vp )*dVdt;

}

//amd
void get_drho_dt(struct planet * pl , double r , double phi , double rho ,  double * drho_dt_sink){

   double r_p = pl->r;
   double p_p = pl->phi;
   double m_p = pl->M;

   // Distance of the cell to the planet:
   double dx = r*cos(phi)-r_p*cos(p_p);
   double dy = r*sin(phi)-r_p*sin(p_p);
   double script_r = sqrt(dx*dx+dy*dy);

   // Assumes r_sink = smoothing length:
   //double eps = pl->eps;  
   //double r_sink = r_sink_factor*eps;
   double r_sink = pl->eps;

   // **Assumes alpha viscosity
   // Must be changed for nu
   double sink_nu = alpha * pow(m_p,0.5) * pow(r_sink,0.5) / pow(mach,2);
 
    // Set accretion timescale to viscous time at rsink
   double t_visc = 2./3. * r_sink * r_sink / sink_nu;
   double t_sink = t_sink_factor * t_visc;

   //*drho_dt_sink = 0.0;

   // If r < r_sink, then drho_dt source term is calculated
   //if (script_r < r_sink){
   //   *drho_dt_sink = rho / t_sink;
   //}
   
   *drho_dt_sink = rho / t_sink * exp(- pow( script_r / r_sink , 4.)); 

}

void planet_sink(struct planet * pl , double * prim , double * cons , double * xp , double * xm , double dVdt ){

   // This is cell data
   double rp = xp[0];
   double rm = xm[0];
   // Radius of the cell
   double r = 0.5*(rp+rm);
   double rho = prim[RHO];
   double pres  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];

   // Phi of cell
   double dphi = get_dp(xp[1],xm[1]);
   double phi = xm[1] + 0.5*dphi;

   double drho_dt_sink;

   // Call get_drho_dt with planet data and cell r, phi
   // rho, pres, viscosity

   get_drho_dt( pl , r , phi , rho , &drho_dt_sink );

   // Update all the conservative variables since they have factors of rho
   // better to use primitive variables on the right hand side, since those 
   // aren't being updated throughout a time step   

   cons[DDD] -= drho_dt_sink*dVdt;
   cons[SRR] -= drho_dt_sink * vr * dVdt;
   cons[TAU] -= .5 * (vr*vr + vp*vp + vz*vz ) * drho_dt_sink * dVdt;
   cons[LLL] -= r * drho_dt_sink * vp * dVdt;
   cons[SZZ] -= drho_dt_sink * vz * dVdt;


}


void planet_RK_copy( struct planet * pl ){
   pl->RK_r     = pl->r;
   pl->RK_phi   = pl->phi;
   pl->RK_M     = pl->M;
   pl->RK_omega = pl->omega;
   pl->RK_vr    = pl->vr;
}

void planet_RK_adjust( struct planet * pl , double RK ){
   pl->r     = (1.-RK)*pl->r     + RK*pl->RK_r;
   pl->phi   = (1.-RK)*pl->phi   + RK*pl->RK_phi;
   pl->M     = (1.-RK)*pl->M     + RK*pl->RK_M;
   pl->omega = (1.-RK)*pl->omega + RK*pl->RK_omega;
   pl->vr    = (1.-RK)*pl->vr    + RK*pl->RK_vr;
}
