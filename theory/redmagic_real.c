double w_clustering_HOD_rm_tomo(double theta, int nz);
double w_clustering_HOD_rm(double theta);
double w_clustering_HOD_rm_max(double theta);//max(w_2h,w_1h) halo exclusion model

double xi_gg (double R, double a);
double xi_gg_1h (double R, double a);
double xi_nl (double R, double a);
/************** real space routines *************/
void w_via_hankel_HOD_rm_tomo(double **xi, double *logthetamin, double *logthetamax, int nz)
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
    //     lntmin   = log(limits.xi_via_hankel_theta_min);
    //     lntmax   = log(limits.xi_via_hankel_theta_max);
    //     dlntheta = (lntmax-lntmin)/(1.0*Ntable.N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = Ntable.N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_cl_HOD_rm_tomo(l,nz);
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

double w_clustering_HOD_rm_tomo(double theta, int nz)
{
  static double MO = 0.;
  static int NZ = -1;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  if (MO != redm.hod[0] || NZ != nz){
    MO = redm.hod[0]; NZ = nz;
    if (table!=0) free_double_matrix(table, 0, 1, 0, Ntable.N_thetaH-1);
    table   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    w_via_hankel_HOD_rm_tomo(table, &logthetamin, &logthetamax, nz);
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH-1.0);
  }
  res = interpol(table[0], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);
  
  return res;
}


/******************** no tomography *********************/
double C_cl_HOD_rm (double k){
  double da,ap,res, array[2];
  array[0] = k; array[1] = 0.;
  return int_gsl_integrate_low_precision(int_for_C_cl_HOD_rm,(void*)array,1./(1+tomo.clustering_zmax[0]),1./(1+tomo.clustering_zmin[0]),NULL,1000);
}

void w_via_hankel_HOD_rm(double **xi, double *logthetamin, double *logthetamax)
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
    //     lntmin   = log(limits.xi_via_hankel_theta_min);
    //     lntmax   = log(limits.xi_via_hankel_theta_max);
    //     dlntheta = (lntmax-lntmin)/(1.0*Ntable.N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = Ntable.N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_cl_HOD_rm(l);
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

double w_clustering_HOD_RM(double theta)
{
  static double MO = 0.;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  if (MO != redm.hod[0]){
    MO = redm.hod[0];
    if (table!=0) free_double_matrix(table, 0, 1, 0, Ntable.N_thetaH-1);
    table   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    w_via_hankel_HOD_rm(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH-1.0);
  }
  res = interpol(table[0], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);
  
  return res;
}
double P_gg_rm_1h(double k,double a)
{
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.,amin =0., amax = 0.;
  static int N_a = 10, N_k_nlin = 50;
  static double M0 = 0.0;
  
  static double **table_P_gg=0;
  
  double klog,val,aa,kk;
  int i,j;
  
  if (M0 != redm.hod[0]){ //extend this by halo model parameters if these become part of the model
    M0 = redm.hod[0];
    if (table_P_gg == 0){
      table_P_gg = create_double_matrix(0, N_a-1, 0, N_k_nlin-1);
    }
    amax = fmin(1./(1+redshift.clustering_zdistrpar_zmin-0.01),0.999);
    amin = 1./(1+redshift.clustering_zdistrpar_zmax+0.01);
    da = (amax - amin)/(N_a-1.);
    logkmin = log(0.1*cosmology.coverH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(N_k_nlin-1.);
    aa= amin;
    for (i=0; i<N_a; i++, aa +=da) {
      klog  = logkmin;
      for (j=0; j<N_k_nlin; j++, klog += dk) {
        kk = exp(klog);
        table_P_gg[i][j]= log(G02_rm(kk,aa)/pow(ngal_rm(aa),2.0));
      }
    }
  }
  klog = log(k);
  val = interpol2d(table_P_gg, N_a, amin, amax, da, a, N_k_nlin, logkmin, logkmax, dk, klog, 0.0, 0.0);
  return exp(val);
}

double int_for_C_cl_HOD_rm_1h(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell       = ar[0]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
    
  res= pf_photoz(1./a-1.,0)*pf_photoz(1./a-1.,0)/(a*a)*hoverh0(a)/fK/fK;
  if (res !=0){res= res*P_gg_rm_1h(k,a);}
  return res;
}

double C_cl_HOD_rm_1h (double k){
  double da,ap,res, array[2];
  array[0] = k;
  return int_gsl_integrate_low_precision(int_for_C_cl_HOD_rm_1h,(void*)array,1./(1+tomo.clustering_zmax[0]),1./(1+tomo.clustering_zmin[0]),NULL,1000);
}

void w_via_hankel_HOD_rm_1h(double **xi, double *logthetamin, double *logthetamax)
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
    //     lntmin   = log(limits.xi_via_hankel_theta_min);
    //     lntmax   = log(limits.xi_via_hankel_theta_max);
    //     dlntheta = (lntmax-lntmin)/(1.0*Ntable.N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = Ntable.N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_cl_HOD_rm_1h(l);
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

double w_clustering_HOD_rm_1h(double theta)
{
  static double MO = 0.;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  if (MO != redm.hod[0]){
    MO = redm.hod[0];
    if (table!=0) free_double_matrix(table, 0, 1, 0, Ntable.N_thetaH-1);
    table   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    w_via_hankel_HOD_rm_1h(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH-1.0);
  }
  res = interpol(table[0], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);
  return res;
}
double int2_for_wreal(double z2,void *params)
{
  double *ar = (double *) params;
  double z1 = ar[1], theta = ar[0];
  double chi1 = chi(1./(1+z1)), chi2 = chi(1./(1+z2));
  double r,xn,br;
  r = sqrt(pow(chi1-chi2,2.0)+pow(0.5*(chi1+chi2)*theta,2.0));
  xn = xi_nl(r,1./(1+0.5*(z1+z2)));
  br = pow(1+1.17*xn,1.49)/pow(1.0+0.69*xn,2.09);
  return pf_photoz(z1,0)*pf_photoz(z2,0)*br*xi_gg(r,1./(1+0.5*(z1+z2)));//*bgal_rm(1./(1+z1))*bgal_rm(1./(1+z2));
}

double int_for_wreal(double z1,void *params)
{
  double *array = (double *) params;
  array[1] = z1;
  return int_gsl_integrate_low_precision(int2_for_wreal,(void*)array,tomo.clustering_zmin[0],tomo.clustering_zmax[0],NULL,1000);
}
double w_real(double theta){
  static double **table = 0;
  static double dlogtheta, logthetamin, logthetamax;
  static int Ntheta_w =200;
  if (table == 0){
    double lgt,array[2];
    int i;
    table   = create_double_matrix(0, 1, 0, Ntheta_w-1);
    logthetamin = log(0.09*constants.arcmin);
    logthetamax = log(120*constants.arcmin);
    dlogtheta = (logthetamax-logthetamin)/(Ntheta_w-1);
    lgt = logthetamin;
    for (i = 0; i < Ntheta_w; i++){
      array[0] = exp(lgt);
      table[0][i] = int_gsl_integrate_low_precision(int_for_wreal,(void*)array,tomo.clustering_zmin[0],tomo.clustering_zmax[0],NULL,1000);
      lgt += dlogtheta;
    }
  }
  return  interpol(table[0], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0)*pow(bgal_rm(1./(1+0.5*tomo.clustering_zmax[0]+0.5*tomo.clustering_zmin[0])),2.0);
}
double w_clustering_HOD_rm(double theta){
  return w_clustering_HOD_RM(theta);
}
double w_clustering_HOD_rm_max(double theta){
  double w1h = w_clustering_HOD_rm_1h(theta);
  return fmax(w_clustering_HOD_RM(theta)-w1h, w1h);
}

/********************* halo exclusion *********************/
/******* halo exclusion routines for YS ************/
double int_for_xi_nl (double k, void *params){
  double *ar = (double *) params;
  return Pdelta(k,ar[1])*k*k*gsl_sf_sinc(k*ar[0]/M_PI);
}

double xi_nl (double R, double a){ // 3D matter correlation function (tabulated at z = 0.5*(tomo.clustering_zmin[0]+tomo.clustering_zmax[0]); scaled with D(z)^2)
  static cosmopara C;
  static double *table;
  static double am = 1.;
  static double logrmin = .0,logrmax=.0,dr=.0;
  if (recompute_cosmo3D(C)){
    printf("tabulating xi_nl\n");
    am = 1./(1+0.5*(tomo.clustering_zmin[0]+tomo.clustering_zmax[0]));
    double d,logr,res, array[2];
    int i;
    long int n;
    
    if (table==0){
      table = create_double_vector(0,Ntable.N_r_3d-1);
      logrmin = log(limits.xi_3d_rmin);
      logrmax = log(limits.xi_3d_rmax);
      dr = (logrmax - logrmin)/(Ntable.N_r_3d-1);
    }
    
    
    logr  = logrmin;
    for (i = 0; i< Ntable.N_r_3d; i++,logr+=dr){
      array[0]= exp(logr);
      array[1] = am;
      n = ceil(limits.k_min_cH0*array[0]/M_PI);
      res = int_gsl_integrate_high_precision(int_for_xi_nl,(void*)array,limits.k_min_cH0,(n*1.0)*M_PI/array[0],NULL,1000);
      d = res;
      while(fabs(d)> 1.e-5*fabs(res) && (n+10.0)*M_PI/array[0] < limits.k_max_cH0){
        d =int_gsl_integrate_medium_precision(int_for_xi_nl,(void*)array,(1.0*n)*M_PI/array[0],(n+10.)*M_PI/array[0],NULL,1000);
        n = n+10;
        res = res+d;
      }
      table[i] = res/(2.0*M_PI*M_PI);
    }
    printf("xi_3D tabulated\n");
    update_cosmopara(&C);
  }
  return interpol(table, Ntable.N_r_3d, logrmin, logrmax,dr, log(R), 1.0, 0.0)*pow(growfac(a)/growfac(am),2.0);
}

double int_P_gg_2h_excl(double lgM, void *para){
  double *array = (double *) para;
  double u, m = exp(lgM),k = array[0],a = array[1];
  array[3] = m;
  u = u_g_rm(k,m,a);
  return m*massfunc(m,a)*(u*n_s_rm(m)+n_c_rm(m))*B1(m,a);
}

double P_gg_2h_excl_old(double k,double r, double a){
  static cosmopara C;
  static double *table, R = -1.;
  static double logkmin = .0,logkmax=.0,dk=.0;
  
  double am = 1./(1+0.5*(tomo.clustering_zmin[0]+tomo.clustering_zmax[0]));
  double M_max_r = m_Delta(r,am);
  if (M_max_r > 1.e+15){return Pdelta(k,a);}
  if (M_max_r < 1.e+11){return 0.;}
  if (recompute_cosmo3D(C) ||R!= r ){
    double logk, array[4]= {0,am,0.,0.}, n2,b2,bb;
    double lgM1, lgM2, M1, M2, kk,res,res1,dm = 0.1*log(10.);
    int i, pm;
    pm = redm.parameterization;
    redm.parameterization = 2;
    logkmin = log(0.5*cosmology.coverH0);
    logkmax = log(20.0*cosmology.coverH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_3d-1);
    
    if (table!=0) free_double_vector(table, 0, Ntable.N_k_3d-1);
    table   = create_double_vector(0,Ntable.N_k_3d-1);
    array[0] = am;
    logk  = logkmin;
    array[0] = am;
    //    b2 = pow(int_gsl_integrate_medium_precision(int_bgal_rm, (void*)array, log(10.)*(11.),log(10.)*(15.),NULL, 5000)*int_gsl_integrate_medium_precision(int_ngal_rm, (void*)array, log(10.)*(11.),log(M_max_r),NULL, 5000)/int_gsl_integrate_medium_precision(int_ngal_rm, (void*)array, log(10.)*(11.),log(10.)*(15.),NULL, 5000),2.0);
    b2 = pow(int_gsl_integrate_medium_precision(int_bgal_rm, (void*)array, log(10.)*(11.),log(M_max_r),NULL, 5000)/int_gsl_integrate_medium_precision(int_ngal_rm, (void*)array, log(10.)*(11.),log(10.)*(15.),NULL, 5000),2.0);
    for (i = 0; i< Ntable.N_k_3d; i++,logk+=dk){
      array[0] = exp(logk);
      kk = exp(logk);
      res = 0;
      //b2 = 0;
      for (lgM1 = log(10.)*11; lgM1 < log(M_max_r-1.e+11); lgM1 +=dm){
        M1 = exp(lgM1);
        res1 = 0.;
        for (lgM2 = log(10.)*11.; lgM2 < log(fmin(M_max_r, m_Delta(radius(M1)-r,am))); lgM2 +=dm){
          M2 = exp(lgM2);
          res1 += M2*massfunc(M2,am)*(u_g_rm(kk,M2,am)*n_s_rm(M2)+n_c_rm(M2))*B1(M2,am)*dm;
        }
        res +=res1*M1*massfunc(M1,am)*(u_g_rm(kk,M1,am)*n_s_rm(M1)+n_c_rm(M1))*B1(M1,am)*dm;
      }
      table[i] = fmin(1.,res/b2);
      //    table[i] = sqrt(fmin(1.0,pow(int_gsl_integrate_low_precision(int_P_gg_2h_excl, (void*)array,log(10.)*(11.),log(M_max_r),NULL, 5000),2.0)/(b2)));
      //      printf("%d/%d: %e %e %e\n", i, Ntable.N_k_3d, kk*cosmology.coverH0,table[i], Pdelta(array[0],am));
    }
    R = r;
    update_cosmopara(&C);
    redm.parameterization = pm;
  }
  if (k <= exp(logkmin)){return Pdelta(k,a)*table[0];}
  if (k >= exp(logkmax)){return Pdelta(k,a)*table[Ntable.N_k_3d-1];}
  return Pdelta(k,a)*interpol(table, Ntable.N_k_3d, logkmin, logkmax,dk, log(k), 1.0, 1.0);
}

double int_for_xi_gg (double k, void *params){
  double *ar = (double *) params;
  return P_gg_rm(k,ar[1])*k*k*gsl_sf_sinc(k*ar[0]/M_PI);
}

double int_for_xi_gg_1h (double k, void *params){
  double *ar = (double *) params;
  return P_gg_rm_1h(k,ar[1])*k*k*gsl_sf_sinc(k*ar[0]/M_PI);
}

/*double xi_gg (double R, double a){ // 3D gg correlation function (tabulated at z = 0.5*(tomo.clustering_zmin[0]+tomo.clustering_zmax[0]); scaled with D(z)^2)
  static cosmopara C;
  static double *table;
  static double am = 1.;
  static double logrmin = .0,logrmax=.0,dr=.0;
  
  if (recompute_cosmo3D(C || am ! = a)){
    am = a;
    double d,logr,res, array[2];
    int i;
    long int n;
    
    if (table==0){
      table = create_double_vector(0,Ntable.N_r_3d-1);
      logrmin = log(limits.xi_3d_rmin);
      logrmax = log(limits.xi_3d_rmax);
      dr = (logrmax - logrmin)/(Ntable.N_r_3d-1);
    }
    
    
    logr  = logrmin;
    for (i = 0; i< Ntable.N_r_3d; i++,logr+=dr){
      array[0]= exp(logr);
      array[1] = am;
      n = 1;
      res = int_gsl_integrate_medium_precision(int_for_xi_gg,(void*)array,1.e-6*M_PI/array[0],(n*1.0)*M_PI/array[0],NULL,1000);
      d = res;
      while(fabs(d)> 1.e-4*fabs(res)){
        d =int_gsl_integrate_medium_precision(int_for_xi_gg,(void*)array,(1.0*n)*M_PI/array[0],(n+2.)*M_PI/array[0],NULL,1000);
        n = n+2;
        res = res+d;
      }
      table[i] = res/(2.0*M_PI*M_PI);
    }
    printf("xi_gg tabulated\n");
    update_cosmopara(&C);
  }
  return interpol(table, Ntable.N_r_3d, logrmin, logrmax,dr, log(R), 1.0, 0.0));
}*/

               
double xi_gg (double R_mpc, double a){ // 3D gg correlation function (tabulated at z = 0.5*(tomo.clustering_zmin[0]+tomo.clustering_zmax[0]); scaled with D(z)^2)
    double d,res, array[2];
    long int n;
    array[0]= R_mpc/cosmology.coverH0;
    array[1] = a;
    n = ceil(limits.k_min_cH0*array[0]/M_PI);
    res = int_gsl_integrate_high_precision(int_for_xi_gg,(void*)array,limits.k_min_cH0,(n*1.0)*M_PI/array[0],NULL,1000);
    d = res;
    while(fabs(d)> 1.e-5*fabs(res) && (n+10.0)*M_PI/array[0] < limits.k_max_cH0){
        d =int_gsl_integrate_medium_precision(int_for_xi_gg,(void*)array,(1.0*n)*M_PI/array[0],(n+10.)*M_PI/array[0],NULL,1000);
        n = n+10;
        res = res+d;
    }
    return res/(2.0*M_PI*M_PI);
}
double xi_gg_1h (double R_mpc, double a){ // 3D gg correlation function (tabulated at z = 0.5*(tomo.clustering_zmin[0]+tomo.clustering_zmax[0]); scaled with D(z)^2)
    double d,res, array[2];
    long int n;
    array[0]= R_mpc/cosmology.coverH0;
    array[1] = a;
    n = ceil(limits.k_min_cH0*array[0]/M_PI);
    res = int_gsl_integrate_high_precision(int_for_xi_gg_1h,(void*)array,limits.k_min_cH0,(n*1.0)*M_PI/array[0],NULL,1000);
    d = res;
    while(fabs(d)> 1.e-5*fabs(res) && (n+10.0)*M_PI/array[0] < limits.k_max_cH0){
        d =int_gsl_integrate_medium_precision(int_for_xi_gg_1h,(void*)array,(1.0*n)*M_PI/array[0],(n+10.)*M_PI/array[0],NULL,1000);
        n = n+10;
        res = res+d;
    }
    return res/(2.0*M_PI*M_PI);
}

double xi_nl2 (double R, double a){ // 3D matter correlation function (tabulated at z = 0.5*(tomo.clustering_zmin[0]+tomo.clustering_zmax[0]); scaled with D(z)^2)
  static cosmopara C;
  static double *table;
  static double am = 1.;
  static double logrmin = .0,logrmax=.0,dr=.0;
  if (recompute_cosmo3D(C)){
    printf("tabulating xi_nl\n");
    am = 1./(1+0.5*(tomo.clustering_zmin[0]+tomo.clustering_zmax[0]));
    double d,logr,res, array[2];
    int i;
    long int n;
    
    if (table==0){
      table = create_double_vector(0,Ntable.N_r_3d-1);
      logrmin = log(limits.xi_3d_rmin);
      logrmax = log(limits.xi_3d_rmax);
      dr = (logrmax - logrmin)/(Ntable.N_r_3d-1);
    }
    
    
    logr  = logrmin;
    for (i = 0; i< Ntable.N_r_3d; i++,logr+=dr){
      array[0]= exp(logr);
      array[1] = am;
      n = ceil(limits.k_min_cH0*array[0]/M_PI);
      res = int_gsl_integrate_high_precision(int_for_xi_nl,(void*)array,limits.k_min_cH0,(n*1.0)*M_PI/array[0],NULL,1000);
      d = res;
      while(fabs(d)> 1.e-5*fabs(res) && (n+10.0)*M_PI/array[0] < limits.k_max_cH0){
        d =int_gsl_integrate_medium_precision(int_for_xi_nl,(void*)array,(1.0*n)*M_PI/array[0],(n+10.)*M_PI/array[0],NULL,1000);
        n = n+10;
        res = res+d;
      }
      table[i] = res/(2.0*M_PI*M_PI);
    }
    printf("xi_3D tabulated\n");
    update_cosmopara(&C);
  }
  return interpol(table, Ntable.N_r_3d, logrmin, logrmax,dr, log(R), 1.0, 0.0)*pow(growfac(a)/growfac(am),2.0);
}

