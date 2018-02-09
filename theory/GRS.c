double f_a(double a);
double Tsqr_EH_no_wiggle(double khoverMPC);
double P_DW(double k, double mu, double a);
double P_obs(double k_ref, double mu_ref, int nz);
double P_obs_mt(double k_ref, double mu_ref, int nz, GRSpara_mt  G, int i, int j);
typedef struct{
   int N_z;
   int N_k;
   int N_mu;
   double k_star; //in h/Mpc
   double k_min; // in h/Mpc
   double k_max; //in h/Mpc
   double f_sky; 
   double z[10];
   double V_z[10]; // in (Mpc/h)^3
   double H_ref[10];
   double DA_ref[10];
   double* datav;
   double* var;
   double* k;
   double* mu;
} GRSpara;
GRSpara GRS;

typedef struct{
  double n_g[10]; // in (h/Mpc)^3
	double b_g[10];
	double sigma_z[10]; // fractional accuracy
	double sigma_p[10]; // in km/s
	double P_shot[10]; // in (Mpc/h)^3
} GRS_galaxy_para;
GRS_galaxy_para GRS_gal;

double logG(double loga,void * params){
     return log(growfac(exp(loga)));
}
double f_a(double a){
  static cosmopara C;
  static double *table_fa;
  static double da = .0, logamin = 1.0,logamax = 1.0;
  double alog,error;
  int i;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);

    if (table_fa==0){
      table_fa  = create_double_vector(0, Ntable.N_a-1);
      logamin = log(limits.a_min);
      logamax = 0.0;
      da = (logamax - logamin)/(Ntable.N_a-1.);
    }
    gsl_function F;
    double result, abserr;
  
    F.function = &logG;
    F.params = 0;
    alog = logamin;
    for (i=0; i<Ntable.N_a; i++, alog += da) {
      gsl_deriv_central (&F,alog, da/100., &result, &abserr);
      table_fa[i]=result;
    }
  }
  return interpol(table_fa, Ntable.N_a, logamin, logamax, da,log(a), 1.0,1.0);

}
double Tsqr_EH_no_wiggle(double khoverMPC)
{
     double q, theta, ommh2, a, s, gamma, L0, C0;
     double tmp;
     double omegam, ombh2, hubble;

     /* other input parameters */
     hubble = cosmology.h0;

     omegam = cosmology.Omega_m;
     ombh2 = cosmology.omb* cosmology.h0 *cosmology.h0 ;

     if(cosmology.omb == 0)
          ombh2 = 0.04 * cosmology.h0 *cosmology.h0 ;

     theta = 2.728 / 2.7;
     ommh2 = omegam * hubble * hubble;
     s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
     a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
          + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
     gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * khoverMPC * s)));
     gamma *= omegam * hubble;
     q = khoverMPC * theta * theta / gamma; 
     L0 = log(2. * exp(1.) + 1.8 * q);
     C0 = 14.2 + 731. / (1. + 62.5 * q);
     tmp = L0 / (L0 + C0 * q * q);
     return (tmp*tmp);
}     

double P_DW(double k, double mu, double a){
  double G2 = pow(growfac(a)/growfac(1.0),2.0);
  double g_mu = G2*(1-mu*mu+mu*mu*pow(1+f_a(a),2.));
  double damping_los = exp(-g_mu*k*k/(GRS.k_star*GRS.k_star));
  double Tsqr = Tsqr_EH_no_wiggle(k)*(1.-damping_los)+Tsqr_EH_wiggle(k)*damping_los;
  return G2*cosmology.sigma_8*cosmology.sigma_8/sigma_r_sqr()*Tsqr*pow(k,cosmology.n_spec);
}
double P_g(double k, double mu, int nz){
  double a =1./(1.+GRS.z[nz]);
 return pow(GRS_gal.b_g[nz]*(1+f_a(a)/GRS_gal.b_g[nz]*mu*mu),2.)*P_DW(k,mu,a);	
}
double P_obs(double k_ref, double mu_ref, int nz){
  double mu, k, k_perp, k_los;
  double a =1./(1.+GRS.z[nz]);

  double H_Href = cosmology.h0*hoverh0(a)/GRS.H_ref[nz];
  double DA_DAref = f_K(a)/cosmology.h0/GRS.DA_ref[nz];
  //calculate (k,mu) corresponding to (k_ref,mu_ref) in current cosmology
  k_perp = k_ref*sqrt(1-mu_ref*mu_ref)/DA_DAref;
  k_los = k_ref*mu_ref*H_Href;
  k = pow(k_los*k_los+k_perp*k_perp,0.5);
  mu = k_los/k;

  //distance dispersion corresponding to the physical velocity dispersion sigma_p[nz]
  double sigma_r_p = (GRS_gal.sigma_p[nz]/2.997e+5)/(cosmology.h0*hoverh0(a)*a);// in c/H0 units
  sigma_r_p *= cosmology.coverH0;
  //distance dispersion corresponding to the redshift dispersion sigma_z[nz]
  double sigma_r_z  = GRS_gal.sigma_z[nz]/(cosmology.h0*hoverh0(a)*a);// in c/H0 units
  sigma_r_z *= cosmology.coverH0;
  //printf("%e %e %e\n",k,sigma_r_z, exp(-pow(k*sigma_r_z,2.)));
  //this expression matches Eq. 4 in Wang et al. 2013
  return pow(DA_DAref,-2.)*H_Href*P_g(k,mu,nz)
  	*exp(-pow(k*mu*sigma_r_z,2.))/(1.+0.5*pow(k*mu*sigma_r_p,2.))
   	+ GRS_gal.P_shot[nz];
}

double P_g_mt(double k, double mu, int nz, GRSpara_mt  G, int i, int j){
  double a =1./(1.+G.z[nz]);
// return G.b_mt[nz][i]*G.b_mt[nz][j]*(1+f_a(a)/G.b_mt[nz][i]*mu*mu)*(1+f_a(a)/G.b_mt[nz][j]*mu*mu)*Pdelta(k*cosmology.coverH0,a)*pow(cosmology.coverH0,3.);  
 return G.b_mt[nz][i]*G.b_mt[nz][j]*(1+f_a(a)/G.b_mt[nz][i]*mu*mu)*(1+f_a(a)/G.b_mt[nz][j]*mu*mu)*p_lin(k*cosmology.coverH0,a)*pow(cosmology.coverH0,3.);  
}
double P_obs_mt(double k_ref, double mu_ref, int nz, GRSpara_mt  G, int i, int j){
//  return P_g_mt(k_ref,mu_ref,nz,G,i,j);
//    *exp(-pow(k*mu*sigma_r_z,2.))/(1.+0.5*pow(k*mu*sigma_r_p,2.));
  double mu, k, k_perp, k_los;
  double a =1./(1.+G.z[nz]);

  double H_Href = cosmology.h0*hoverh0(a)/G.H_ref[nz];
  double DA_DAref = f_K(a)/cosmology.h0/G.DA_ref[nz];
  //calculate (k,mu) corresponding to (k_ref,mu_ref) in current cosmology
  k_perp = k_ref*sqrt(1-mu_ref*mu_ref)/DA_DAref;
  k_los = k_ref*mu_ref*H_Href;
  k = pow(k_los*k_los+k_perp*k_perp,0.5);
  mu = k_los/k;

  //distance dispersion corresponding to the physical velocity dispersion sigma_p[nz]
  double sigma_r_p = (G.sigma_p/2.997e+5)/(cosmology.h0*hoverh0(a)*a);// in c/H0 units
  sigma_r_p *= cosmology.coverH0;
  //distance dispersion corresponding to the redshift dispersion sigma_z[nz]
  double sigma_r_z  = G.sigma_z/(cosmology.h0*hoverh0(a)*a);// in c/H0 units
  sigma_r_z *= cosmology.coverH0;

  //this expression matches Eq. 4 in Wang et al. 2013
  return pow(DA_DAref,-2.)*H_Href*P_g_mt(k,mu,nz,G,i,j)
    *exp(-pow(k*mu*sigma_r_z,2.))/(1.+0.5*pow(k*mu*sigma_r_p,2.));
}


