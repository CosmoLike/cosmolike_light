/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
double cov_NG_shear_shear_real(double theta1, double theta2, int z1,int z2, int z3, int z4, int pm1, int pm2); //pm1 = pm2 = 1 for C++, pm1=pm2 = 0 for C--, pm1 = 1, pm2 = 0 for C+-
double cov_G_shear_shear_real(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4, int pm1, int pm2);//pm1 = pm2 = 1 for C++, pm1=pm2 = 0 for C--, pm1 = 1, pm2 = 0 for C+-
double cov_G_shear_rebin(double thetamin_i, double thetamax_i,double thetamin_j,double thetamax_j, int z1,int z2,int z3,int z4,int pm1,int pm2); //Version of Gaussian cov calculation for wide bins
double cov_NG_shear_rebin(double thetamin_i, double thetamax_i,double thetamin_j,double thetamax_j, int z1,int z2,int z3,int z4,int pm1,int pm2);

double cov_NG_cl_cl_real(double theta1, double theta2, int z1,int z2, int z3, int z4);
double cov_G_cl_cl_real(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4);
double cov_G_cl_cl_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j, double thetamax_j, int z1,int z2, int z3, int z4);

double cov_NG_gl_gl_real(double theta1, double theta2, int z1,int z2, int z3, int z4);
double cov_G_gl_gl_real(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4);
double cov_G_gl_gl_real_rebin(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1l,int z1s, int z2l, int z2s);

double cov_NG_cl_shear_real(double theta1, double theta2, int z1,int z2, int z3, int z4, int pm); //z1,z2 clustering bins; z3,z4 shear bins
double cov_G_cl_shear_real(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4, int pm);
double cov_G_cl_shear_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j,double thetamax_j, int z1,int z2, int z3, int z4, int pm);

double cov_NG_cl_gl_real(double theta1, double theta2, int z1,int z2, int zl, int zs);//z1,z2 clustering bins; zl,zs g-g lensing bins
double cov_G_cl_gl_real(double theta1, double theta2, double Dtheta, int z1,int z2, int zl, int zs);
double cov_G_cl_gl_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j,double thetamax_j, int z1,int z2, int zl, int zs);

double cov_NG_gl_shear_real(double theta1, double theta2, int zl,int zs, int z3, int z4, int pm);
double cov_G_gl_shear_real(double theta1, double theta2, double Dtheta, int zl,int zs, int z3, int z4,int pm);
double cov_G_gl_shear_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j,double thetamax_j,int zl,int zs, int z3, int z4,int pm);


/************************* covariance routines for angular correlation functions *********************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

double C_gl_tomo_all(double l, int ni, int nj)  //slower version of G-G lensing power spectrum, lens bin ni, source bin nj - tabulated for all lens-source combinations without overlap criterion
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1 = 0.;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo_all(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
    if (table==0){
      if (tomo.shear_Nbin> 10){printf("tomo.shear_Nbin too large for look-up table in C_gl_tomo_all\nEXIT\n");exit(1);}
      table   = create_double_matrix(0, tomo.clustering_Nbin*10-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    int i,j,k;
    double llog;
    
    for (k=0; k<tomo.clustering_Nbin; k++) {
      for (j=0; j<tomo.shear_Nbin; j++) {
        llog = logsmin;
        for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          table[k*10+j][i]= log(C_gl_tomo_nointerp(exp(llog),k,j));
        }
      }
    }
    
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
  }
  f1 = exp(interpol(table[10*ni+nj], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.,1.));
  if (isnan(f1)){f1 = 0;}
  return f1;
  // return C_gl_tomo(l,ni,nj);
}
/**************** shear-shear *********************/
double bin_cov_NG_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 50;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(2.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_shear_shear_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,1.0,1.0));}
  return res;
}

double int2_for_cov_NG_shear(double l2,void *params){
  double *ar = (double *) params;
  double l1,tri = 0.,res =0;
  int n1,n2,n3,n4,j0;
  l1= ar[0];
  j0= (int)ar[8];
  
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4];
  //printf("hsv %le %le\n",l1,l2);
  //hsv = HSV_shear_shear_tomo(l1,l2,n1,n2,n3,n4)*survey.area*survey.area_conversion_factor;
  tri= bin_cov_NG_shear_shear_tomo(l1,l2,n1,n2,n3,n4);
  //t1h = project_tri_1h_cov_shear_shear_tomo(l1,l2,n1,n2,n3,n4);
  //t2h = project_tri_2h_cov_shear_shear_tomo(l1,l2,n1,n2,n3,n4);
  if (j0==1) res=(tri)*l2*gsl_sf_bessel_J0(l2*ar[6]);
  if (j0==0) res=(tri)*l2*gsl_sf_bessel_Jnu(4.,l2*ar[6]);
  return res;
}

double int_for_cov_NG_shear(double l1, void *params){
  double *array = (double *) params,res=0.,result = 1.,res_fin = 0.;
  int j0;
  array[0] = l1;
  j0= (int)array[8];
  unsigned int n =1;
  double x2 =0, x1 = 20.; //x1 = l_min for NG Covariance integration
  while (x2 <= x1){ //find first root of J0/4(l*theta2) with l > 20
    if (j0==1) x2 = gsl_sf_bessel_zero_J0 (n)/array[6];
    if (j0==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/array[6];
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4, strongly shot noise dominated afterwards
    result=int_gsl_integrate_low_precision(int2_for_cov_NG_shear,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    if (j0==1) x2 = gsl_sf_bessel_zero_J0 (n)/array[6];
    if (j0==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/array[6];
    n+=2;
  }
  if ((int)array[7]==1) res_fin=res*l1*gsl_sf_bessel_J0(l1*array[5]);
  if ((int)array[7]==0) res_fin=res*l1*gsl_sf_bessel_Jnu(4.,l1*array[5]);
  return res_fin;
}


double int_for_cov_G_shear(double l, void *params){
  double *ar = (double *) params;
  double JJ,C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = (int) ar[1];
  n2 = (int) ar[2];
  n3 = (int) ar[3];
  n4 = (int) ar[4];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_shear_tomo(l,n1,n3);
  C24 = C_shear_tomo(l,n2,n4);
  C14 = C_shear_tomo(l,n1,n4);
  C23 = C_shear_tomo(l,n2,n3);
  
  if (n1 == n3){N13= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  
  if ((int)ar[7] == 1){JJ = gsl_sf_bessel_J0(l*ar[5]);}
  else{JJ = gsl_sf_bessel_Jnu(4.,l*ar[5]);}
  
  if ((int)ar[8] == 1){JJ *= gsl_sf_bessel_J0(l*ar[6]);}
  else{JJ *= gsl_sf_bessel_Jnu(4.,l*ar[6]);}
  
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*l*JJ;
}

double int_for_cov_G_shear_nonoise(double l, void *params){ //no mixed term!
  double *ar = (double *) params;
  double JJ,C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = (int) ar[1];
  n2 = (int) ar[2];
  n3 = (int) ar[3];
  n4 = (int) ar[4];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_shear_tomo(l,n1,n3);
  C24 = C_shear_tomo(l,n2,n4);
  C14 = C_shear_tomo(l,n1,n4);
  C23 = C_shear_tomo(l,n2,n3);

  if ((int)ar[7] == 1){JJ = gsl_sf_bessel_J0(l*ar[5]);}
  else{JJ = gsl_sf_bessel_Jnu(4.,l*ar[5]);}

  if ((int)ar[8] == 1){JJ *= gsl_sf_bessel_J0(l*ar[6]);}
  else{JJ *= gsl_sf_bessel_Jnu(4.,l*ar[6]);}
  
  return (C13*C24+ C14*C23)*l*JJ;
}

double cov_NG_shear_shear_real(double theta1, double theta2, int z1,int z2, int z3, int z4, int pm1, int pm2){ //pm1 = pm2 = 1 for C++, pm1=pm2 = 0 for C--, pm1 = 1, pm2 = 0 for C+-
  double array[9],res = 0, result = 1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3; array[4] = (double) z4;
  array[5] = theta1; array[6] = theta2; array[7] = (double) pm1;array[8] = (double) pm2;
  
  unsigned int n=1;
  double x2 =0, x1 = 20.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0/4(l*theta1) with l > 20
    if (pm1==1) x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    if (pm1==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/theta1;
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) &&x1<1.e+4){ //integrate up to l_max = 1.e+4, strongly shot noise dominated afterwards
    result=int_gsl_integrate_low_precision(int_for_cov_NG_shear,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    if (pm1==1) x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    if (pm1==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/theta1;
    n+=2;
  }
  return res/(4.0*M_PI*M_PI);
}


//Gaussian part of the shear covariance Cov(xi_pm1(theta1,z1,z2),xi_pm2(theta2,z3,z4)); theta1, theta2, Dtheta in radian, z1..z4 number of tomography bins
//Only one Dtheta is needed as it enters only in the (shot noise)^2 term, which requires theta1 = theta2,
double cov_G_shear_shear_real(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4, int pm1, int pm2){
  //pm1 = pm2 = 1 for C++, pm1=pm2 = 0 for C--, pm1 = 1, pm2 = 0 for C+-
  //printf("%le %le\n",theta1,theta2);
  double N =0., array[9],res = 0., result =1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;  array[7] = (double) pm1;array[8] = (double) pm2;
  if (z1 ==z3 && z2 ==z4 && fabs(theta1-theta2)< 0.1*Dtheta && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    N = pow(survey.sigma_e,4.0)/(2.*M_PI*theta1*Dtheta*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
  }
  if (z1 ==z4 && z2 ==z3 && fabs(theta1-theta2)< 0.1*Dtheta && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    N += pow(survey.sigma_e,4.0)/(2*M_PI*theta1*Dtheta*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
  }
  array[5] = theta1; array[6] = theta2;
  unsigned int n = 1;
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0/4(l*theta1) with l > 1
    if (pm1==1) x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    if (pm1==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/theta1;
    //printf("%le %d\n",x2,n);
    
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<5.e+4){
  //integrate up to l_max = 1.e6
    //printf("%le %le\n",x1,x2);
    result=int_gsl_integrate_medium_precision(int_for_cov_G_shear,(void*)array, x1, x2 ,NULL,1000);
    res = res+result;
    x1 = x2;
    if (pm1==1) x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    if (pm1==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/theta1;
    n+=2;
  }
//  if (N)  printf("%e %e\n",res/(2.0*M_PI*survey.area*survey.area_conversion_factor),2.*N);

  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor)+2.*N;
}

/****************** clustering ***************************/
double bin_cov_NG_cl_cl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 50;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(2.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_cl_cl_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
  }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,1.0,1.0));}
  return res;
}
double bin_cov_NG_cl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 50;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(2.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_cl_shear_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,1.0,1.0));}
  return res;
}

double bin_cov_NG_cl_gl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 50;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(2.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_cl_gl_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,1.0,1.0));}
  return res;
}

double int2_for_cov_NG_cl_cl(double l2,void *params){
  double *ar = (double *) params;
  double l1,tri =0;
  int n1,n2,n3,n4;
  l1= ar[0];
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4];
  tri = bin_cov_NG_cl_cl_tomo(l1,l2,n1,n2,n3,n4);
  return tri*l2*gsl_sf_bessel_J0(l2*ar[6]);
}
double int2_for_cov_NG_cl_shear(double l2,void *params){
  double *ar = (double *) params;
  double l1,tri =0;
  int n1,n2,n3,n4,j0;
  l1= ar[0];
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4];j0 = (int)ar[7];
  tri = bin_cov_NG_cl_shear_tomo(l1,l2,n1,n2,n3,n4);
  if (j0 ==1){return tri*l2*gsl_sf_bessel_J0(l2*ar[6]);}
  return tri*l2*gsl_sf_bessel_Jnu(4.,l2*ar[6]);
}
double int2_for_cov_NG_cl_gl(double l2,void *params){
  double *ar = (double *) params;
  double l1,tri =0;
  int n1,n2,n3,n4,pm;
  l1= ar[0];
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4];
  tri = bin_cov_NG_cl_gl_tomo(l1,l2,n1,n2,n3,n4);
  return tri*l2*gsl_sf_bessel_Jnu(2.,l2*ar[6]);
}

double int_for_cov_NG_cl_cl(double l1, void *params){
  double *array = (double *) params,res =0.,result = 1;
  array[0] = l1;
  unsigned int n =1;
  double x2 =0, x1 = 20.; //x1 = l_min for NG Covariance integration
  while (x2 <= x1){ //find first root of J0(l*theta2) with l > 20
    x2 = gsl_sf_bessel_zero_J0 (n)/array[6];
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4, strongly shot noise dominated afterwards
    result=int_gsl_integrate_low_precision(int2_for_cov_NG_cl_cl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_J0 (n)/array[6];
    n+=2;
  }
  return res*l1*gsl_sf_bessel_J0(l1*array[5]);
}

double int_for_cov_NG_cl_shear(double l1, void *params){
  double *array = (double *) params,res=0.,result = 1.,res_fin = 0.;
  int j0;
  array[0] = l1;
  j0= (int)array[7];
  unsigned int n =1;
  double x2 =0, x1 = 20.; //x1 = l_min for NG Covariance integration
  while (x2 <= x1){ //find first root of J0/4(l*theta2) with l > 20
    if (j0==1) x2 = gsl_sf_bessel_zero_J0 (n)/array[6];
    if (j0==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/array[6];
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4, strongly shot noise dominated afterwards
    result=int_gsl_integrate_low_precision(int2_for_cov_NG_cl_shear,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    if (j0==1) x2 = gsl_sf_bessel_zero_J0 (n)/array[6];
    if (j0==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/array[6];
    n+=2;
  }
  return res*l1*gsl_sf_bessel_J0(l1*array[5]);
}

double int_for_cov_NG_cl_gl(double l1, void *params){
  double *array = (double *) params,res =0.,result = 1.;
  array[0] = l1;
  unsigned int n=1;
  double x2 =0, x1 = 20.; //x1 = l_min for NG Covariance integration
  while (x2 <= x1){ //find first root of J2(l*theta2) with l > 20
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/array[6];
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4, strongly shot noise dominated at higher l
    result=int_gsl_integrate_low_precision(int2_for_cov_NG_cl_gl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/array[6];
    n+=2;
  }
  return res*l1*gsl_sf_bessel_J0(l1*array[5]);
}

//Non-Gaussian part of the clustering covariance Cov(w(theta1,z1,z2),w(theta2,z3,z4)); theta1, theta2 in radian, z1..z4 number of tomography bins

double cov_NG_cl_cl_real(double theta1, double theta2, int z1,int z2, int z3, int z4){
  double array[7],res = 0., result = 1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;
  array[5] = theta1; array[6] = theta2;
  
  unsigned int n=1;
  double x2 =0, x1 = 20.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0(l*theta1) with l > 20
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4
    result=int_gsl_integrate_low_precision(int_for_cov_NG_cl_cl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n+=2;
  }
  return res/(4.0*M_PI*M_PI);
}
double cov_NG_cl_shear_real(double theta1, double theta2, int z1,int z2, int z3, int z4, int pm){
  double array[8],res = 0., result = 1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;
  array[5] = theta1; array[6] = theta2; array[7] = (double)pm;
  
  unsigned int n=1;
  double x2 =0, x1 = 20.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0(l*theta1) with l > 20
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4
    result=int_gsl_integrate_low_precision(int_for_cov_NG_cl_shear,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n+=2;
  }
  return res/(4.0*M_PI*M_PI);
}
double cov_NG_cl_gl_real(double theta1, double theta2, int z1,int z2, int z3, int z4){
  double array[7],res = 0., result = 1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;
  array[5] = theta1; array[6] = theta2;
  
  unsigned int n=1;
  double x2 =0, x1 = 20.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0(l*theta1) with l > 20
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4
    result=int_gsl_integrate_low_precision(int_for_cov_NG_cl_gl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n+=2;
  }
  return res/(4.0*M_PI*M_PI);
}



//Gaussian part of the clustering covariance Cov(w(theta1,z1,z2),w(theta2,z3,z4)); theta1, theta2, Dtheta in radian, z1..z4 number of tomography bins
//Only one Detheta is needed as it enters only in the (shot noise)^2 term, which requires theta1 = theta2,
double int_for_cov_G_cl_cl(double l, void *params){
  double *ar = (double *) params;
  double C13, C14, C23, C24, N13 =0., N14=0., N23=0., N24=0.;
  int n1,n2,n3,n4;
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4];
  C13 = C_cl_tomo(l,ar[1],ar[3]);C24 = C_cl_tomo(l,ar[2],ar[4]);
  C14 = C_cl_tomo(l,ar[1],ar[4]);C23 = C_cl_tomo(l,ar[2],ar[3]);
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= 1./(nlens(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= 1./(nlens(n2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*l*gsl_sf_bessel_J0(l*ar[5])*gsl_sf_bessel_J0(l*ar[6]);
}

double int_for_cov_G_cl_shear(double l, void *params){
  double *ar = (double *) params;
  double C13, C14, C23, C24, N13 =0., N14=0., N23=0., N24=0.;
  int n1,n2,n3,n4, j0;
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4]; j0= (int)ar[7];
  C13 = C_gl_tomo_all(l,n1,n3);C24 = C_gl_tomo_all(l,n2,n4);
  C14 = C_gl_tomo_all(l,n1,n4);C23 = C_gl_tomo_all(l,n2,n3);
  
  if (j0 ==1)return (C13*C24+  C14*C23)*l*gsl_sf_bessel_J0(l*ar[5])*gsl_sf_bessel_J0(l*ar[6]);
  if (j0 ==0)return (C13*C24+  C14*C23)*l*gsl_sf_bessel_J0(l*ar[5])*gsl_sf_bessel_Jnu(4.,l*ar[6]);
  printf("int_for_cov_G_cl_shear: pm = %d invalid value\n(use j0 =1 for xi_+, j0 = 0 for xi_i)\n", j0);
  exit(0);
}

double int_for_cov_G_cl_gl(double l, void *params){
  double *ar = (double *) params;
  double C13, C14, C23, C24, N13 =0., N14=0., N23=0., N24=0.;
  int n1,n2,n3,n4;
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4];
  C13 = C_cl_tomo(l,n1,n3);C24 = C_gl_tomo_all(l,n2,n4);
  C14 = C_gl_tomo_all(l,n1,n4);C23 = C_cl_tomo(l,n2,n3);
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*l*gsl_sf_bessel_J0(l*ar[5])*gsl_sf_bessel_Jnu(2.,l*ar[6]);
}


double cov_G_cl_cl_real(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4){
  double N =0, array[7],res = 0, result = 1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;
  if (z1 ==z3 && z2 ==z4 && fabs(theta1-theta2)< 0.1*Dtheta){
    N = 1./(2.*M_PI*theta1*Dtheta*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
  }
  if (z1 ==z4 && z2 ==z3 && fabs(theta1-theta2)< 0.1*Dtheta){
    N += 1./(2.*M_PI*theta1*Dtheta*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
  }
  array[5] = theta1; array[6] = theta2;
  double y = theta1;
  unsigned int n = 1;
  if (theta2 > y) {y = theta2;} //sample l integral according to the stronger oscillating Bessel function;integrate piece wise from root n to root n+2
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0(l*min(theta1,theta2)) with l > 10
    x2 = gsl_sf_bessel_zero_J0 (n)/y;
    n++;
  }
  while (fabs(result) > 1.e-5*fabs(res) && x1<5.e+4){ //integrate up to l_max = 5.e+4
    result=int_gsl_integrate_medium_precision(int_for_cov_G_cl_cl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_J0 (n)/y;
    n+=2;
  }
  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor)+N;
}
double cov_G_cl_cl_real_nonoise(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4){
  double array[7],res = 0, result = 1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;
  array[5] = theta1; array[6] = theta2;
  double y = theta1;
  unsigned int n = 1;
  if (theta2 > y) {y = theta2;} //sample l integral according to the stronger oscillating Bessel function;integrate piece wise from root n to root n+2
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0(l*min(theta1,theta2)) with l > 10
    x2 = gsl_sf_bessel_zero_J0 (n)/y;
    n++;
  }
  while (fabs(result) > 1.e-5*fabs(res) && x1<5.e+4){ //integrate up to l_max = 5.e+4
    result=int_gsl_integrate_medium_precision(int_for_cov_G_cl_cl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_J0 (n)/y;
    n+=2;
  }
  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor);
}

double cov_G_cl_shear_real(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4, int pm){
  double N =0, array[8],res = 0, result = 1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;array[7]= (double) pm;
  array[5] = theta1; array[6] = theta2;
  unsigned int n = 1;
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0(l*theta1) with l > 10
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n++;
  }
  while (fabs(result) > 1.e-4*fabs(res) && x1<5.e+4){ //integrate up to l_max = 5.e+4
    result=int_gsl_integrate_medium_precision(int_for_cov_G_cl_shear,(void*)array, x1, x2 ,NULL,1000);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n+=20;
  }
  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor);
}

double cov_G_cl_gl_real(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4){
  double N =0, array[7],res = 0, result = 1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;
  array[5] = theta1; array[6] = theta2;
  unsigned int n = 1;
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0(l*theta1) with l > 10
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n++;
  }
  while (fabs(result) > 1.e-5*fabs(res) && x1<5.e+4){ //integrate up to l_max = 5.e+4
    result=int_gsl_integrate_medium_precision(int_for_cov_G_cl_gl,(void*)array, x1, x2 ,NULL,1000);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    n+=2;
  }
  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor);
}


/****************** galaxy-galaxy lensing ***************************/
double bin_cov_NG_gl_gl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 50;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(2.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_gl_gl_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,1.0,1.0));}
  return res;
}
double bin_cov_NG_gl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 50;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(2.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_gl_shear_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,1.0,1.0));}
  return res;
}

double int2_for_cov_NG_gl_gl(double l2,void *params){
  double *ar = (double *) params;
  double l1,tri;
  int n1,n2,n3,n4;
  l1= ar[0];
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4];
  tri = bin_cov_NG_gl_gl_tomo(l1,l2,n1,n2,n3,n4);
  return tri*l2*gsl_sf_bessel_Jn(2,l2*ar[6]);
}
double int2_for_cov_NG_gl_shear(double l2,void *params){
  double *ar = (double *) params;
  double l1,tri =0;
  int n1,n2,n3,n4,j0;
  l1= ar[0];
  n1 = (int) ar[1];n2 = (int) ar[2];n3 = (int) ar[3];n4 = (int) ar[4];j0 = (int)ar[7];
  tri = bin_cov_NG_gl_shear_tomo(l1,l2,n1,n2,n3,n4);
  if (j0 ==1){return tri*l2*gsl_sf_bessel_J0(l2*ar[6]);}
  return tri*l2*gsl_sf_bessel_Jnu(4.,l2*ar[6]);
}

double int_for_cov_NG_gl_gl(double l1, void *params){
  double *array = (double *) params,res =0.,result = 1.;
  array[0] = l1;
  unsigned int n=1;
  double x2 =0, x1 = 20.; //x1 = l_min for NG Covariance integration
  while (x2 <= x1){ //find first root of J2(l*theta2) with l > 20
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/array[6];
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4, strongly shot noise dominated at higher l
    result=int_gsl_integrate_low_precision(int2_for_cov_NG_gl_gl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/array[6];
    n+=2;
  }
  return res*l1*gsl_sf_bessel_Jn(2.,l1*array[5]);
}
double int_for_cov_NG_gl_shear(double l1, void *params){
  double *array = (double *) params,res=0.,result = 1.,res_fin = 0.;
  int j0;
  array[0] = l1;
  j0= (int)array[7];
  unsigned int n =1;
  double x2 =0, x1 = 20.; //x1 = l_min for NG Covariance integration
  while (x2 <= x1){ //find first root of J0/4(l*theta2) with l > 20
    if (j0==1) x2 = gsl_sf_bessel_zero_J0 (n)/array[6];
    if (j0==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/array[6];
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4, strongly shot noise dominated afterwards
    result=int_gsl_integrate_low_precision(int2_for_cov_NG_gl_shear,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    if (j0==1) x2 = gsl_sf_bessel_zero_J0 (n)/array[6];
    if (j0==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/array[6];
    n+=2;
  }
  return res*l1*gsl_sf_bessel_Jn(2,l1*array[5]);
}

double int_for_cov_G_gl_gl(double l, void *params){
  double *ar = (double *) params;
  double C13, C14, C23, C24, N13 =0., N14=0.,N24=0., N23=0.;
  C13 = C_cl_tomo(l,(int) ar[1],(int) ar[3]);
  C24 = C_shear_tomo(l,(int) ar[2],(int) ar[4]);
  C14 = C_gl_tomo_all(l,(int) ar[1],(int) ar[4]);
  C23 = C_gl_tomo_all(l,(int) ar[3],(int) ar[2]);
  if ((int) ar[1] == (int) ar[3]){N13= 1./(nlens((int) ar[1])*survey.n_gal_conversion_factor);}
  if ((int) ar[2] == (int) ar[4]){N24= pow(survey.sigma_e,2.0)/(2.0*nsource((int) ar[2])*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*l*gsl_sf_bessel_Jn(2,l*ar[5])*gsl_sf_bessel_Jn(2,l*ar[6]);
}

double int_for_cov_G_gl_shear(double l, void *params){
  double *ar = (double *) params;
  int j0 = (int) ar[7];
  double C13, C14, C23, C24, N23 =0., N24=0.,N13=0.,N14=0.;
  C13 = C_gl_tomo_all(l,(int) ar[1],(int) ar[3]);
  C24 = C_shear_tomo(l,(int) ar[2],(int) ar[4]);
  C14 = C_gl_tomo_all(l,(int) ar[1],(int) ar[4]);
  C23 = C_shear_tomo(l,(int) ar[2],(int) ar[3]);
  if ((int) ar[2] == (int) ar[4]){N24= pow(survey.sigma_e,2.0)/(2.0*nsource((int) ar[2])*survey.n_gal_conversion_factor);}
  if ((int) ar[2] == (int) ar[3]){N23= pow(survey.sigma_e,2.0)/(2.0*nsource((int) ar[2])*survey.n_gal_conversion_factor);}
  
  if (j0 ==1){return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*l*gsl_sf_bessel_Jn(2,l*ar[5])*gsl_sf_bessel_J0(l*ar[6]);}
  if (j0 ==0){return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*l*gsl_sf_bessel_Jn(2,l*ar[5])*gsl_sf_bessel_Jn(4,l*ar[6]);}
  printf("int_for_cov_G_gl_shear: pm = %d invalid value\n(use j0 =1 for xi_+, j0 = 0 for xi_i)\n", j0);
  exit(0);
}


//Non-Gaussian part of the G-G lensing covariance Cov(gamma_t(theta1,z1l,z1s),gamma_t(theta2,z2l,z2s)); theta1, theta2, Dtheta in radian, z1l,z2l: lens redshift bins, z1s,z2s: source redshift bins
double cov_NG_gl_gl_real(double theta1, double theta2, int z1l,int z1s, int z2l, int z2s){
  double array[7],res = 0., result = 1.;
  array[1] = (double) z1l; array[2] = (double) z1s;array[3] = (double) z2l;array[4] = (double) z2s;
  array[5] = theta1; array[6] = theta2;
  unsigned int n = 1;
  double x2 =0, x1 = 20.; //x1 = l_min for NG Covariance integration
  while (x2 <= x1){ //find first root of J2(l*theta1) with l > 20
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/theta1;
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<5.e+4){ //integrate up to l_max = 1.e+4, dominated by shot noise at higher l..
    result=int_gsl_integrate_low_precision(int_for_cov_NG_gl_gl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/theta1;
    n+=2;
  }
  return res/(4.0*M_PI*M_PI);
}
double cov_NG_gl_shear_real(double theta1, double theta2, int z1l,int z1s, int z3, int z4, int pm){
  double array[8],res = 0., result = 1.;
  array[1] = (double) z1l; array[2] = (double) z1s;array[3] = (double) z3;array[4] = (double) z4; array[7]=(double) pm;
  array[5] = theta1; array[6] = theta2;
  unsigned int n = 1;
  double x2 =0, x1 = 20.; //x1 = l_min for NG Covariance integration
  while (x2 <= x1){ //find first root of J2(l*theta1) with l > 20
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/theta1;
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<1.e+4){ //integrate up to l_max = 1.e+4, dominated by shot noise at higher l..
    result=int_gsl_integrate_low_precision(int_for_cov_NG_gl_shear,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/theta1;
    n+=2;
  }
  return res/(4.0*M_PI*M_PI);
}

//Gaussian part of the G-G lensing covariance Cov(gamma_t(theta1,z1l,z1s),gamma_t(theta2,z2l,z2s)); theta1, theta2, Dtheta in radian, z1l,z2l: lens redshift bins, z1s,z2s: source redshift bins
//Only one Dtheta is needed as it enters only in the (shot noise)^2 term, which requires theta1 = theta2
double cov_G_gl_gl_real(double theta1, double theta2, double Dtheta, int z1l,int z1s, int z2l, int z2s){
  double N=0., array[7],res = 0., result = 1.;
  array[1] = (double) z1l; array[2] = (double) z1s;array[3] = (double) z2l;array[4] = (double) z2s;
  array[5] = theta1; array[6] = theta2;
  if (z1l ==z2l && z1s ==z2s && fabs(theta1-theta2)< 0.1*Dtheta){
    N = pow(survey.sigma_e,2.0)/(2.0*2.*M_PI*theta1*Dtheta*nlens(z1l)*nsource(z2s)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
  }
  double y = theta1;
  unsigned int n = 1;
  if (theta2 > y) {y = theta2;} //sample l integral according to the stronger oscillating Bessel function;integrate piece wise from root n to root n+2
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J2(l*min(theta1,theta2)) with l > 10
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/y;
    n++;
    
  }
  while (fabs(result) > 1.e-5*fabs(res) && x1<5.e+4){ //integrate up to l_max = 2.e+4
    result=int_gsl_integrate_medium_precision(int_for_cov_G_gl_gl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/y;
    n+=2;
  }
  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor)+N;
}

double cov_G_ggl_no_shot_noise(double theta1, double theta2, double Dtheta, int z1l,int z1s, int z2l, int z2s){
  double N=0., array[7],res = 0., result = 1.;
  array[1] = (double) z1l; array[2] = (double) z1s;array[3] = (double) z2l;array[4] = (double) z2s;
  array[5] = theta1; array[6] = theta2;
  double y = theta1;
  unsigned int n = 1;
  if (theta2 > y) {y = theta2;} //sample l integral according to the stronger oscillating Bessel function;integrate piece wise from root n to root n+2
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J2(l*min(theta1,theta2)) with l > 10
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/y;
    n++;
    
  }
  while (fabs(result) > 1.e-5*fabs(res) && x1<5.e+4){ //integrate up to l_max = 2.e+4
    result=int_gsl_integrate_medium_precision(int_for_cov_G_gl_gl,(void*)array, x1, x2 ,NULL,512);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/y;
    n+=2;
  }
  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor);
}



double cov_G_gl_shear_real(double theta1, double theta2, double Dtheta, int z1l,int z1s, int z3, int z4, int pm){
  double N=0., array[8],res = 0., result = 1.;
  array[1] = (double) z1l; array[2] = (double) z1s;array[3] = (double) z3;array[4] = (double) z4;array[7]=(double) pm;
  array[5] = theta1; array[6] = theta2;
  unsigned int n = 1;
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J2(l*min(theta1,theta2)) with l > 10
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/theta1;
    n++;
    
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<5.e+4){ //integrate up to l_max = 2.e+4
    result=int_gsl_integrate_medium_precision(int_for_cov_G_gl_shear,(void*)array, x1, x2 ,NULL,1000);
    res = res+result;
    x1 = x2;
    x2 = gsl_sf_bessel_zero_Jnu (2.,n)/theta1;
    n+=100;
  }
  
  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor);
}








/*********** rebin routines for shear shear **********/
double cov_G_shear_no_shot_noise(double theta1, double theta2, double Dtheta, int z1,int z2, int z3, int z4, int pm1, int pm2){
 double N =0., array[9],res = 0., result =1.;
  array[1] = (double) z1; array[2] = (double) z2;array[3] = (double) z3;array[4] = (double) z4;  array[7] = (double) pm1;array[8] = (double) pm2;
  array[5] = theta1; array[6] = theta2;
  unsigned int n = 1;
  double x2 =0, x1 = 1.; //x1 = l_min for Covariance integration, OK for typical angular scales of interest
  while (x2 <= x1){ //find first root of J0/4(l*theta1) with l > 1
    if (pm1==1) x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    if (pm1==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/theta1;
    //printf("%le %d\n",x2,n);
    
    n++;
  }
  while (fabs(result) > 1.e-3*fabs(res) && x1<5.e+4){
  //integrate up to l_max = 1.e6
    //printf("%le %le\n",x1,x2);
    result=int_gsl_integrate_medium_precision(int_for_cov_G_shear,(void*)array, x1, x2 ,NULL,1000);
    res = res+result;
    x1 = x2;
    if (pm1==1) x2 = gsl_sf_bessel_zero_J0 (n)/theta1;
    if (pm1==0) x2 = gsl_sf_bessel_zero_Jnu (4.,n)/theta1;
    n+=2;
  }
  return res/(2.0*M_PI*survey.area*survey.area_conversion_factor);
}
double cov_NG_shear_rebin(double thetamin_i, double thetamax_i,double thetamin_j,double thetamax_j, int z1,int z2,int z3,int z4,int pm1,int pm2){
  int Nsub_G = 3;
  int ii,jj;
  double ti,tj,dti, dtj;
  double cG = 0.;
  dti = (thetamax_i-thetamin_i)/(double)Nsub_G;
  dtj = (thetamax_j-thetamin_j)/(double)Nsub_G;
  for (ii = 0; ii < Nsub_G; ii++){
    ti = 2./3.*(pow(thetamin_i+(ii+1.)*dti,3.)-pow(thetamin_i+(ii+0.)*dti,3.))/(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));
    for (jj = 0; jj < Nsub_G; jj++){
      tj = 2./3.*(pow(thetamin_j+(jj+1.)*dtj,3.)-pow(thetamin_j+(jj+0.)*dtj,3.))/(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
      cG+= cov_NG_shear_shear_real(ti,tj, z1,z2,z3,z4,pm1,pm2)*(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.))*(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
    }
  }
  return cG/((pow(thetamax_i,2.)-pow(thetamin_i,2.))*(pow(thetamax_j,2.)-pow(thetamin_j,2.)));
}


double cov_G_shear_shotnoise(double thetamin_i, double thetamax_i,double thetamin_j,double thetamax_j, int z1,int z2,int z3,int z4,int pm1,int pm2){
  double N= 0.;
  if (z1 ==z3 && z2 ==z4 && fabs(thetamax_i-thetamax_j)< 0.1*(thetamax_j-thetamin_j) && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    N = pow(survey.sigma_e,4.0)/(M_PI*(pow(thetamax_i,2.)-pow(thetamin_i,2.0))*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
  }
  if (z1 ==z4 && z2 ==z3 && fabs(thetamax_i-thetamax_j)< 0.1*(thetamax_j-thetamin_j) && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    N += pow(survey.sigma_e,4.0)/(M_PI*(pow(thetamax_i,2.)-pow(thetamin_i,2.0))*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
  }
  return N;
}



double cov_G_cl_shear_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j, double thetamax_j, int z1,int z2, int z3, int z4, int pm){
  double cG = 0.;
  int Nsub_G =3;
  double dti = (thetamax_i-thetamin_i)/(double)Nsub_G;
  double dtj = (thetamax_j-thetamin_j)/(double)Nsub_G;
  for (int ii = 0; ii < Nsub_G; ii++){
    double ti = 2./3.*(pow(thetamin_i+(ii+1.)*dti,3.)-pow(thetamin_i+(ii+0.)*dti,3.))/(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));
    for (int jj = 0; jj < Nsub_G; jj++){
      double tj = 2./3.*(pow(thetamin_j+(jj+1.)*dtj,3.)-pow(thetamin_j+(jj+0.)*dtj,3.))/(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
      cG+= cov_G_cl_shear_real(ti,tj, dti,z1,z2,z3,z4,pm)*(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.))*(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
    }
  }
  return cG/((pow(thetamax_i,2.)-pow(thetamin_i,2.))*(pow(thetamax_j,2.)-pow(thetamin_j,2.)));
}

double cov_G_gl_shear_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j, double thetamax_j, int zl,int zs, int z3, int z4, int pm){
  double cG = 0.;
  int Nsub_G =2;
  double dti = (thetamax_i-thetamin_i)/(double)Nsub_G;
  double dtj = (thetamax_j-thetamin_j)/(double)Nsub_G;
  for (int ii = 0; ii < Nsub_G; ii++){
    double ti = 2./3.*(pow(thetamin_i+(ii+1.)*dti,3.)-pow(thetamin_i+(ii+0.)*dti,3.))/(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));
    for (int jj = 0; jj < Nsub_G; jj++){
      double tj = 2./3.*(pow(thetamin_j+(jj+1.)*dtj,3.)-pow(thetamin_j+(jj+0.)*dtj,3.))/(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
      cG+= cov_G_gl_shear_real(ti,tj, dti,zl,zs,z3,z4,pm)*(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.))*(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
    }
  }
  return cG/((pow(thetamax_i,2.)-pow(thetamin_i,2.))*(pow(thetamax_j,2.)-pow(thetamin_j,2.)));
}

double cov_G_cl_gl_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j, double thetamax_j, int z1,int z2, int z3, int z4){
  double cG = 0.;
  int Nsub_G =3;
  double dti = (thetamax_i-thetamin_i)/(double)Nsub_G;
  double dtj = (thetamax_j-thetamin_j)/(double)Nsub_G;
  for (int ii = 0; ii < Nsub_G; ii++){
    double ti = 2./3.*(pow(thetamin_i+(ii+1.)*dti,3.)-pow(thetamin_i+(ii+0.)*dti,3.))/(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));
    for (int jj = 0; jj < Nsub_G; jj++){
      double tj = 2./3.*(pow(thetamin_j+(jj+1.)*dtj,3.)-pow(thetamin_j+(jj+0.)*dtj,3.))/(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
      cG+= cov_G_cl_gl_real(ti,tj, dti,z1,z2,z3,z4)*(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.))*(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
    }
  }
  return cG/((pow(thetamax_i,2.)-pow(thetamin_i,2.))*(pow(thetamax_j,2.)-pow(thetamin_j,2.)));
}

double cov_G_cl_cl_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j, double thetamax_j, int z1,int z2, int z3, int z4){
  double N=0., array[9],cG = 0., result = 1.;
  int Nsub_G =3;
  if (z1 ==z3 && z2 ==z4 && fabs(thetamin_i-thetamin_j)<1.e-7){
    N = 1./(M_PI*(pow(thetamax_i,2.)-pow(thetamin_i,2.0))*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
    Nsub_G = 5;
  }
  if (z1 ==z4 && z2 ==z3 && fabs(thetamin_i-thetamin_j)<1.e-7){
    N += 1./(M_PI*(pow(thetamax_i,2.)-pow(thetamin_i,2.0))*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    Nsub_G = 5;
  }
  double dti = (thetamax_i-thetamin_i)/(double)Nsub_G;
  double dtj = (thetamax_j-thetamin_j)/(double)Nsub_G;
  for (int ii = 0; ii < Nsub_G; ii++){
    double ti = 2./3.*(pow(thetamin_i+(ii+1.)*dti,3.)-pow(thetamin_i+(ii+0.)*dti,3.))/(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));
    for (int jj = 0; jj < Nsub_G; jj++){
      double tj = 2./3.*(pow(thetamin_j+(jj+1.)*dtj,3.)-pow(thetamin_j+(jj+0.)*dtj,3.))/(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
      cG+= cov_G_cl_cl_real_nonoise(ti,tj, dti,z1,z2,z3,z4)*(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.))*(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
    }
  }
  return cG/((pow(thetamax_i,2.)-pow(thetamin_i,2.))*(pow(thetamax_j,2.)-pow(thetamin_j,2.)))+N;
}

double cov_G_gl_gl_real_rebin(double thetamin_i, double thetamax_i,double thetamin_j, double thetamax_j, int z1l,int z1s, int z2l, int z2s){
  double N=0., array[9],cG = 0., result = 1.;
  int Nsub_G =3;
  if (z1l ==z2l && z1s ==z2s && fabs(thetamin_i-thetamin_j)<1.e-7){
    N = pow(survey.sigma_e,2.0)/(2.0*M_PI*(thetamax_i*thetamax_i-thetamin_i*thetamin_i)*nlens(z1l)*nsource(z2s)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    Nsub_G = 5;
  }
  double dti = (thetamax_i-thetamin_i)/(double)Nsub_G;
  double dtj = (thetamax_j-thetamin_j)/(double)Nsub_G;
  for (int ii = 0; ii < Nsub_G; ii++){
    double ti = 2./3.*(pow(thetamin_i+(ii+1.)*dti,3.)-pow(thetamin_i+(ii+0.)*dti,3.))/(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));
    for (int jj = 0; jj < Nsub_G; jj++){
      double tj = 2./3.*(pow(thetamin_j+(jj+1.)*dtj,3.)-pow(thetamin_j+(jj+0.)*dtj,3.))/(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
      cG+= cov_G_ggl_no_shot_noise(ti,tj, dti,z1l,z1s,z2l,z2s)*(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.))*(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
    }
  }
  return cG/((pow(thetamax_i,2.)-pow(thetamin_i,2.))*(pow(thetamax_j,2.)-pow(thetamin_j,2.)))+N;
}
double cov_G_shear_rebin(double thetamin_i, double thetamax_i,double thetamin_j,double thetamax_j, int z1,int z2,int z3,int z4,int pm1,int pm2){ //now without shotnoise!
  int Nsub_G = 3;
  int ii,jj;
  double ti,tj,dti, dtj;
  double cG = 0., N= 0.;
  if (z1 ==z3 && z2 ==z4 && fabs(thetamax_i-thetamax_j)< 0.1*(thetamax_j-thetamin_j) && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    N = pow(survey.sigma_e,4.0)/(M_PI*(pow(thetamax_i,2.)-pow(thetamin_i,2.0))*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
    Nsub_G = 5;
  }
  if (z1 ==z4 && z2 ==z3 && fabs(thetamax_i-thetamax_j)< 0.1*(thetamax_j-thetamin_j) && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    N += pow(survey.sigma_e,4.0)/(M_PI*(pow(thetamax_i,2.)-pow(thetamin_i,2.0))*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    Nsub_G = 5;
  }
  
  dti = (thetamax_i-thetamin_i)/(double)Nsub_G;
  dtj = (thetamax_j-thetamin_j)/(double)Nsub_G;
  for (ii = 0; ii < Nsub_G; ii++){
    ti = 2./3.*(pow(thetamin_i+(ii+1.)*dti,3.)-pow(thetamin_i+(ii+0.)*dti,3.))/(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));
    for (jj = 0; jj < Nsub_G; jj++){
      tj = 2./3.*(pow(thetamin_j+(jj+1.)*dtj,3.)-pow(thetamin_j+(jj+0.)*dtj,3.))/(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
      cG+= cov_G_shear_no_shot_noise(ti,tj, dti,z1,z2,z3,z4,pm1,pm2)*(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.))*(pow(thetamin_j+(jj+1.)*dtj,2.)-pow(thetamin_j+(jj+0.)*dtj,2.));
    }
  }
  //if (N)  printf("%e %e\n",cG/((pow(thetamax_i,2.)-pow(thetamin_i,2.))*(pow(thetamax_j,2.)-pow(thetamin_j,2.))),2.*N);
  return cG/((pow(thetamax_i,2.)-pow(thetamin_i,2.))*(pow(thetamax_j,2.)-pow(thetamin_j,2.)))+ 2.0*N;
}