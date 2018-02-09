
/*********** Limber approximation integrands for angular power spectra *************/
double int_for_C_magnification_magnification(double a, void *params);
double int_for_C_shear_magnification(double a, void *params);
double int_for_C_position_magnification(double a, void *params);

/********* angular power spectra **********/
double C_magnification_magnification_tomo(double l, int ni, int nj); //magnification tomography power spectra
double C_magnification_magnification(double l);//magnification non-tomography power spectra
double C_shear_mag_tomo(double l, int ni, int nj); //shear-magnification tomography power spectra

double C_shear_mag(double l); //shear-magnification non-tomography power spectra
double C_pos_mag_tomo(double l, int ni, int nj); //position-magnification tomography power spectra
double C_pos_mag(double l); //position-magnification non-tomography power spectra
/*************** hankel transformation routines for angular correlation functions ***************/
void xi_via_hankel_magnification_magnification_tomo(double **xi, double *logthetamin, double *logthetamax, int ni, int nj);
void xi_via_hankel_shear_magnification_tomo(double **xi, double *logthetamin, double *logthetamax, int ni, int nj);
void xi_via_hankel_position_magnification_tomo(double **xi, double *logthetamin, double *logthetamax, int ni, int nj);

/**************** angular correlation functions ***************/
double xi_magnification_magnification_tomo(double theta, int ni, int nj);
//magnification magnification tomography 2PCF galaxies in bins ni, nj
double xi_shear_magnification_tomo(double theta, int ni, int nj);
//shear magnification tomography 2PCF galaxies in bins ni, nj
double xi_position_magnification_tomo(double theta, int ni, int nj);
//magnification galaxy position tomography 2PCF galaxies in bins ni, nj

/*********** Limber approximation integrands for angular power spectra *************/

double int_for_C_magnification_magnification(double a, void *params)
{
  double *ar = (double *) params;
  double res,hoverh0,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res= g_tomo_magnification(a,(int) ar[0])/(a*a)*g_tomo_magnification(a,(int) ar[1])/(a*a)/hoverh0;
  res= res*Pdelta(k,a); //k in units H0/c
  return res;
}

double int_for_C_shear_magnification(double a, void *params)
{
  double *ar = (double *) params;
  double res,hoverh0,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res= g_tomo(a,(int) ar[0])/(a*a)*g_tomo_magnification(a,(int) ar[1])/(a*a)/hoverh0;
  res= res*Pdelta(k,a); //k in units H0/c
  return res;
}

double int_for_C_position_magnification(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= bgal_z(1./a-1., (int)(ar[0]))*pf_photoz(1./a-1.,(int) (ar[0]))*g_tomo_magnification(a,(int) ar[1])/(a*a*a)/fK;
  res= res*Pdelta(k,a);
  return res;
}


/*********** angular power spectra *************/

double C_magnification_magnification_tomo(double s, int ni, int nj) //shear tomography power spectra
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double N_SPEC  = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double SIGMA_8 = -123.;
  
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog;
  int i,j,k,l1,l2;
  
  double array[3],res = 0.0;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8){
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(Ntable.N_ell - 1.);
    
    if (table!=0) free_double_matrix(table, 0, tomo.magnification_Nbin*tomo.magnification_Nbin, 0, Ntable.N_ell-1);
    table   = create_double_matrix(0, tomo.magnification_Nbin*tomo.magnification_Nbin, 0, Ntable.N_ell-1);
    
    array[0] = -1; array[1] = -1;
    slog = logsmin;
    for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
      array[2]= exp(slog);
      res =int_gsl_integrate_medium_precision(int_for_C_magnification_magnification,(void*)array,1./(redshift.magnification_zdistrpar_zmax+1.),0.999999,NULL,1000);
      table[0][i]= log(9.*SQR(cosmology.Omega_m )*(res));
    }
    
    for (k=0; k<tomo.magnification_Nbin; k++) {
      array[0]=(double) k;
      for (j=k; j<tomo.magnification_Nbin; j++) {
	      array[1]=(double) j;
	      l1 = k*tomo.magnification_Nbin+j+1;
        l2 = j*tomo.magnification_Nbin+k+1;
	      slog = logsmin;
	      for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
	        array[2]= exp(slog);
	        res =int_gsl_integrate_medium_precision(int_for_C_magnification_magnification,(void*)array,1./(redshift.magnification_zdistrpar_zmax+1.),0.999999,NULL,1000);
	        table[l1][i]= log(9.*SQR(cosmology.Omega_m )*(res));
	        table[l2][i]= log(9.*SQR(cosmology.Omega_m )*(res));
        }
      }
    }
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
  }
  slog = log(s);
  if (nj == -1 && ni == -1){
    f1 = exp(interpol(table[0], Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));}
  else { f1 = exp(interpol(table[ni*tomo.magnification_Nbin+nj+1], Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));}
  if (isnan(f1)){f1 = 0;}
  return f1;
}

double C_magnification_magnification(double s) 
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double N_SPEC  = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double SIGMA_8 = -123.;
  
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog;
  int i;
  
  double array[3],res = 0.0;
  array[0]= -1;
  array[1]= -1;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8){
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(Ntable.N_ell - 1.);
    
    if (table!=0) free_double_vector(table, 0, Ntable.N_ell-1);
    table   = create_double_vector(0, Ntable.N_ell-1);
    slog = logsmin;
	for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
	  array[2]= exp(slog);
	  res =int_gsl_integrate_medium_precision(int_for_C_magnification_magnification,(void*)array,1./(redshift.magnification_zdistrpar_zmax+1.),0.999999,NULL,1000);
	  table[i]= log(9.*SQR(cosmology.Omega_m )*(res));
    }
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
  }
  slog = log(s);
  f1 = exp(interpol(table, Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));
  if (isnan(f1)){f1 = 0;}
  return f1;
}



double C_shear_mag_tomo(double s, int ni, int nj) //shear magnification tomography power spectra
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double N_SPEC  = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double SIGMA_8 = -123.;
  
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog;
  int i,j,k,l;
  
  double array[3],res = 0.0;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8){
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(Ntable.N_ell - 1.);
    
    if (table!=0) free_double_matrix(table, 0, tomo.shear_Nbin*tomo.magnification_Nbin, 0, Ntable.N_ell-1);
    table   = create_double_matrix(0, tomo.shear_Nbin*tomo.magnification_Nbin, 0, Ntable.N_ell-1);
    
    array[0] = -1; array[1] = -1;
    slog = logsmin;
    for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
      array[2]= exp(slog);
      res =int_gsl_integrate_medium_precision(int_for_C_shear_magnification,(void*)array,1./(redshift.magnification_zdistrpar_zmax+1.),0.999999,NULL,1000);
      table[0][i]= log(9./2.*SQR(cosmology.Omega_m )*(res));
    }

    for (k=0; k<tomo.shear_Nbin; k++) {
      array[0]=(double) k;
      for (j=0; j<tomo.magnification_Nbin; j++) {
        array[1]=(double) j;
	      l = k*tomo.shear_Nbin+j+1;
	      slog = logsmin;
	      for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
	         array[2]= exp(slog);
	         res =int_gsl_integrate_medium_precision(int_for_C_shear_magnification,(void*)array,1./(redshift.shear_zdistrpar_zmax+1.),0.999999,NULL,1000);
	         table[l][i]= log(9./2.*SQR(cosmology.Omega_m )*(res));
	      }
	}
    }
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
  }
  slog = log(s);
  if (nj == -1 && ni == -1){
    f1 = exp(interpol(table[0], Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));}
  else {f1 = exp(interpol(table[ni*tomo.shear_Nbin+nj+1], Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));}
  if (isnan(f1)){f1 = 0;}
  return f1;
}


double C_shear_mag(double s) 
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double N_SPEC  = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double SIGMA_8 = -123.;
  
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog;
  int i;
  
  double array[3],res = 0.0;
  array[0]= -1;
  array[1]= -1;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8){
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(Ntable.N_ell - 1.);
    
    if (table!=0) free_double_vector(table, 0, Ntable.N_ell-1);
    table   = create_double_vector(0, Ntable.N_ell-1);
    slog = logsmin;
	for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
	  array[2]= exp(slog);
	  res =int_gsl_integrate_medium_precision(int_for_C_shear_magnification,(void*)array,1./(redshift.magnification_zdistrpar_zmax+1.),0.999999,NULL,1000);
	  table[i]= log(9./2.*SQR(cosmology.Omega_m )*(res));
    }
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
  }
  slog = log(s);
  f1 = exp(interpol(table, Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));
  if (isnan(f1)){f1 = 0;}
  return f1;
}


double C_pos_mag_tomo(double s, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double N_SPEC  = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double SIGMA_8 = -123.;
  
  static double BIAS[11]   = {-123.,-123.,-123.,-123.,-123.,-123.,-123.,-123.,-123.,-123.,-123.};
  
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double slog,f1;
  int i,j,k,l;
  double array[3],res = 0.0;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BIAS[ni+1] != gbias.b[ni+1][0] )
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(Ntable.N_ell - 1.);
    slog = logsmin;
    if (table!=0) free_double_matrix(table, 0, tomo.shear_Nbin*tomo.clustering_Nbin, 0, Ntable.N_ell-1);
    table   = create_double_matrix(0, tomo.shear_Nbin*tomo.clustering_Nbin, 0, Ntable.N_ell-1);
  
    array[0] = -1; array[1] = -1;
    for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
	    array[2] = exp(slog);
	    res =int_gsl_integrate_medium_precision(int_for_C_position_magnification,(void*)array,limits.a_min,0.9999,NULL,1000);
	    table[0][i]= log(3.*cosmology.Omega_m*res);
	  }
    BIAS[0]   = gbias.b[0][0];

    for (k=0; k<tomo.clustering_Nbin; k++) {
      array[0]=(double) k;
      for (j=0; j<tomo.magnification_Nbin; j++) {
        l = k*tomo.clustering_Nbin + j;
        array[1]=(double) j;
        slog = logsmin;
        for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
          array[2] = exp(slog);
          res =int_gsl_integrate_medium_precision(int_for_C_position_magnification,(void*)array,limits.a_min,0.9999,NULL,1000);
          table[l][i]= log(3.*cosmology.Omega_m*res);
        }
      }
      BIAS[k+1]   = gbias.b[k+1][0];
    }
    
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    
  }
  slog = log(s);
  if (nj == -1 && nj == -1){f1 = exp(interpol(table[0], Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));}
  else {f1 = exp(interpol(table[ni*tomo.clustering_Nbin + nj + 1], Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));}
  if (isnan(f1)){f1 = 0;}
  return f1;
}


double C_pos_mag(double s) 
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double N_SPEC  = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double SIGMA_8 = -123.;
  
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog;
  int i;
  
  double array[3],res = 0.0;
  array[0]= -1;
  array[1]= -1;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8){
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(Ntable.N_ell - 1.);
    
    if (table!=0) free_double_vector(table, 0, Ntable.N_ell-1);
    table   = create_double_vector(0, Ntable.N_ell-1);
    slog = logsmin;
	for (i=0; i<Ntable.N_ell; i++, slog+=ds) {
	  array[2]= exp(slog);
	  res =int_gsl_integrate_medium_precision(int_for_C_shear_magnification,(void*)array,1./(redshift.magnification_zdistrpar_zmax+1.),0.999999,NULL,1000);
	  table[i]= log(3.*cosmology.Omega_m*(res));
    }
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
  }
  slog = log(s);
  f1 = exp(interpol(table, Ntable.N_ell, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0));
  if (isnan(f1)){f1 = 0;}
  return f1;
}

/*************** hankel transformation routines for angular correlation functions ***************/

void xi_via_hankel_magnification_magnification_tomo(double **xi, double *logthetamin, double *logthetamax, int ni, int nj)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -123.0, loglmin, dlnl, lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i;
  lP   = fftw_malloc(Ntable.N_thetaH*sizeof(double));
  f_lP = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(Ntable.N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(Ntable.N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-123.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*Ntable.N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = Ntable.N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_magnification_magnification_tomo(l, ni, nj);
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
  arg[1] = 0;   /* order of Bessel function */
  /* perform the convolution, negative sign for kernel (complex conj.!) */
  for(i=0; i<Ntable.N_thetaH/2+1; i++) {
    kk = 2*constants.pi*i/(dlnl*Ntable.N_thetaH);
    hankel_kernel_FT(kk, &kernel, arg, 2);
    conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
    conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
  }
  /* force Nyquist- and 0-frequency-components to be double */
  conv[0][1] = 0;
  conv[Ntable.N_thetaH/2][1] = 0;
  /* go back to double space, i labels log-bins in theta */
  fftw_execute(plan1);
  for(i=0; i<Ntable.N_thetaH; i++) {
    t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
    xi[0][Ntable.N_thetaH-i-1] = lP[i]/(t*2*constants.pi*Ntable.N_thetaH);
  }
  
  *logthetamin = (nc-Ntable.N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}



void xi_via_hankel_shear_magnification_tomo(double **xi, double *logthetamin, double *logthetamax, int ni, int nj)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -123.0, loglmin, dlnl, lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i;
  lP   = fftw_malloc(Ntable.N_thetaH*sizeof(double));
  f_lP = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(Ntable.N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(Ntable.N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-123.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*Ntable.N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = Ntable.N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_shear_mag_tomo(l,ni,nj);
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
  arg[1] = 2;   /* order of Bessel function */
  /* perform the convolution, negative sign for kernel (complex conj.!) */
  for(i=0; i<Ntable.N_thetaH/2+1; i++) {
    kk = 2*constants.pi*i/(dlnl*Ntable.N_thetaH);
    hankel_kernel_FT(kk, &kernel, arg, 2);
    conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
    conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
  }
  /* force Nyquist- and 0-frequency-components to be double */
  conv[0][1] = 0;
  conv[Ntable.N_thetaH/2][1] = 0;
  /* go back to double space, i labels log-bins in theta */
  fftw_execute(plan1);
  for(i=0; i<Ntable.N_thetaH; i++) {
    t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
    xi[0][Ntable.N_thetaH-i-1] = lP[i]/(t*2*constants.pi*Ntable.N_thetaH);
  }
  
  *logthetamin = (nc-Ntable.N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}



void xi_via_hankel_position_magnification_tomo(double **xi, double *logthetamin, double *logthetamax, int ni, int nj)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -123.0, loglmin, dlnl, lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i;
  lP   = fftw_malloc(Ntable.N_thetaH*sizeof(double));
  f_lP = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(Ntable.N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(Ntable.N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-123.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*Ntable.N_thetaH-1);
    //dlntheta = (lntmax-lntmin)/(1.0*Ntable.N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = Ntable.N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_pos_mag_tomo(l,ni,nj);
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
  arg[1] = 0;   /* order of Bessel function */
  /* perform the convolution, negative sign for kernel (complex conj.!) */
  for(i=0; i<Ntable.N_thetaH/2+1; i++) {
    kk = 2*constants.pi*i/(dlnl*Ntable.N_thetaH);
    hankel_kernel_FT(kk, &kernel, arg, 2);
    conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
    conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
  }
  /* force Nyquist- and 0-frequency-components to be double */
  conv[0][1] = 0;
  conv[Ntable.N_thetaH/2][1] = 0;
  /* go back to double space, i labels log-bins in theta */
  fftw_execute(plan1);
  for(i=0; i<Ntable.N_thetaH; i++) {
    t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
    xi[0][Ntable.N_thetaH-i-1] = lP[i]/(t*2*constants.pi*Ntable.N_thetaH);
  }
  
  *logthetamin = (nc-Ntable.N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}

/**************** angular correlation functions ***************/
/******************** all angles in radian! *******************/


double xi_magnification_magnification_tomo(double theta, int ni, int nj)
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double SIGMA_8 = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double N_SPEC  = -123.;
  
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double **tab;
  double res;
  
  int i,j,k;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 ) 
  {
    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    OMB =cosmology.omb;
    H0 =cosmology.h0;
    
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    
    if (table!=0) free_double_matrix(table, 0, tomo.magnification_Nbin*tomo.magnification_Nbin-1, 0, Ntable.N_thetaH-1);
    table   = create_double_matrix(0, tomo.magnification_Nbin*tomo.magnification_Nbin-1, 0, Ntable.N_thetaH-1);
    for (i = 0; i < tomo.clustering_Nbin; i++){
      for (j= 0; j < tomo.clustering_Nbin; j++){
	if (table!=0) free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
	tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);     xi_via_hankel_magnification_magnification_tomo(tab, &logthetamin, &logthetamax,ni,nj);
	for (k = 0; k < Ntable.N_thetaH; k++){
	  table[i*tomo.magnification_Nbin+j][k] = tab[0][k];	
	}
	free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
      }
    }
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH-1.0);	
  }
  res = interpol(table[ni*tomo.magnification_Nbin+nj], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res; 
}



double xi_shear_magnification_tomo(double theta, int ni, int nj)
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double SIGMA_8 = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double N_SPEC  = -123.;
  
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double **tab;
  double res;
  
  int i,j,k;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 ) 
  {
    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    OMB =cosmology.omb;
    H0 =cosmology.h0;
    
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    if (table!=0) free_double_matrix(table, 0, tomo.magnification_Nbin*tomo.magnification_Nbin-1, 0, Ntable.N_thetaH-1);
    table   = create_double_matrix(0, tomo.magnification_Nbin*tomo.magnification_Nbin-1, 0, Ntable.N_thetaH-1);
    for (i = 0; i < tomo.shear_Nbin; i++){
      for (j= 0; j < tomo.magnification_Nbin; j++){
	if (table!=0) free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
	tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);      xi_via_hankel_shear_magnification_tomo(tab, &logthetamin, &logthetamax,ni,nj);
	for (k = 0; k < Ntable.N_thetaH; k++){
	  table[i*tomo.shear_Nbin+j][k] = tab[0][k];	
	}
	free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
      }
    }
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH-1.0);	
  }
  res = interpol(table[ni*tomo.magnification_Nbin+nj], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res; 
}



double xi_position_magnification_tomo(double theta, int ni, int nj)
{
  static double OMEGA_M = -123.;
  static double OMEGA_V = -123.;
  static double W0= -123.;
  static double WA = -123.;
  static double N_SPEC  = -123.;
  static double OMB   = -123.;
  static double H0   = -123.;
  static double SIGMA_8 = -123.;
  
  static double BIAS   = -123.;
  static double RCORR   = -123.;
  
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double **tab;
  double res;
  
  int i,j,k;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 )
  {
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    
    if (table!=0) free_double_matrix(table, 0, tomo.magnification_Nbin*tomo.clustering_Nbin-1, 0, Ntable.N_thetaH-1);
    table   = create_double_matrix(0, tomo.magnification_Nbin*tomo.clustering_Nbin-1, 0, Ntable.N_thetaH-1);
    for (i = 0; i < tomo.clustering_Nbin; i++){
      for (j= 0; j < tomo.magnification_Nbin; j++){
	if (table!=0) free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
	tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);  
	xi_via_hankel_position_magnification_tomo(tab, &logthetamin, &logthetamax,ni,nj);
	for (k = 0; k < Ntable.N_thetaH; k++){
	  table[i*tomo.clustering_Nbin+j][k] = tab[0][k];	
	}
	free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
      }
    }
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH-1.0);
  }
  res = interpol(table[ni*tomo.clustering_Nbin+nj], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  
  return res; 
}


////////////// Start magnification routines //////////////////////
////////////// Start magnification routines //////////////////////
/*double g_tomo_magnification(double a, int zbin) // for tomography bin zbin
{
  static cosmopara C;
  static nuisancepara N;
  static double **table;

  static double da = 0.0;
  double aa;
  int i,j;
  double array[2];
  
  if (recompute_expansion(C) ||  recompute_zphot_magnification(N))
    {
    if (table==0)  table = create_double_matrix(0, tomo.magnification_Nbin, 0, Ntable.N_a-1);
    da = (0.999999-1./(redshift.magnification_zdistrpar_zmax+1.))/(Ntable.N_a-1);
    
    for (j=-1;j<tomo.magnification_Nbin;j++) {
      array[0]=(double) j;
      aa = 1./(redshift.magnification_zdistrpar_zmax+1.);
      for (i=0;i<Ntable.N_a;i++,aa+=da) {
        array[1] = aa;
        table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g_tomo,(void*)array,1./(redshift.magnification_zdistrpar_zmax+1.),aa,NULL,4000);
      }
      
    }
    
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (a<1./(redshift.magnification_zdistrpar_zmax+1.)) return 0.0;
  return interpol(table[zbin+1], Ntable.N_a, 1./(redshift.magnification_zdistrpar_zmax+1.), 0.999999, da, a, 1.0, 1.0);
}


double int_for_g_t_magnification(double aprime,void *params)
{
  double chi1, chi_prime,val;
  double *ar = (double *) params;
  int zbin= (int) ar[0];
  chi1 = chi(ar[1]);
  chi_prime = chi(aprime);
  val=zdistr_magnification_photoz(1./aprime-1.,zbin)*f_K(chi_prime-chi1)/f_K(chi_prime)/(aprime*aprime);
  return val;
}

double zdistr_magnification_photoz(double zz,int j) //returns n(ztrue | j), works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
{
  static double **table;
  static double da = 0.0;
  static int tab = 0;
  double static zhisto_max,zhisto_min;
  
  double array[3];
  int i,k;
  if (tab ==0){
    if (redshift.magnification_photoz != 0){
      printf ("magnifition_photoz != 0 currently not supported. Exit.\n");
      exit(1);
    }
    table   = create_double_matrix(0, tomo.magnification_Nbin, 0, redshift.magnification_histogram_zbins-1);
    zhisto_max =redshift.magnification_zdistrpar_zmax;
    zhisto_min = redshift.magnification_zdistrpar_zmin;
    da = (zhisto_max-zhisto_min)/(redshift.magnification_histogram_zbins-1);
    double norm,zi;
    for (i = -1; i< tomo.magnification_Nbin; i++){
      if (i < 0){
        array[0] = redshift.magnification_zdistrpar_zmin;
        array[1] = redshift.magnification_zdistrpar_zmax;
      }
      
      if (i>= 0 && i < tomo.magnification_Nbin){
        array[0] =tomo.magnification_zmin[i];
        array[1] = tomo.magnification_zmax[i];
      }
      norm = int_gsl_integrate_medium_precision(int_for_zdistr_magnification_mock, (void*)array, array[0],array[1],NULL, 1024);
      for (k = 0,zi = zhisto_min;k<redshift.magnification_histogram_zbins; k++,zi+=da){
        table[i+1][k] = int_for_zdistr_magnification_mock(zi,(void*)array)/norm;
      }

    }
    tab = 1;
  }
  if (zz >= zhisto_max || zz < zhisto_min) return 0.0;
  return interpol(table[j+1], redshift.magnification_histogram_zbins, zhisto_min, zhisto_max, da, zz, 1.0, 1.0);
}



double int_for_zdistr_magnification_mock(double z, void *params)
{
  static double *tab;
  char filename[200];
  FILE *ein;
  
  double static zhisto_max,zhisto_min,dz;
  
  if (tab==0){
    double *z_v,space1,space2;
    int i;
    tab=create_double_vector(0, redshift.magnification_histogram_zbins-1);
    z_v=create_double_vector(0, redshift.magnification_histogram_zbins-1);
    sprintf(filename,"%s",redshift.magnification_REDSHIFT_FILE);
    ein=fopen(filename,"r");
    EXIT_MISSING_FILE(ein, "int_for_zdistr_magnification_mock", filename);
    
    if (!ein) {fprintf(stderr,"Could not open magnification_mock file %s\n", filename);exit(1);}
    for (i=0;i<redshift.magnification_histogram_zbins;i++){
      fscanf(ein,"%le %le %le %le\n",&z_v[i],&space1,&space2,&tab[i]);
    }
    fclose(ein);
    dz = z_v[1]-z_v[0];
    zhisto_max=z_v[redshift.magnification_histogram_zbins-1]+dz;
    zhisto_min=z_v[0];
    free_double_vector(z_v,0,redshift.magnification_histogram_zbins-1);
  }
  
  if ((z>=zhisto_min) &&(z<zhisto_max)){
    return tab[(int)floor((z-zhisto_min)/dz)];
  }
  return 0.;
}

*/

/////////////// End magnification routines///////////
/////////////// End magnification routines///////////


