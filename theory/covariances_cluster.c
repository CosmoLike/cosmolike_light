//#cluster covariance: counts + cluster lensing
double cov_G_clg_clg(double l, double delta_l,int nzc1,int nN1, int nzs1, int nzc2, int nN2, int nzs2);
double cov_NG_cgl_clg(double l1,int nzc1,int nN1, int nzs1, double l2, int nzc2, int nN2, int nzs2);
double cov_cgl_N(double l,int nzc1,int nN1, int nzs,int nzc2, int nN2);
double cov_N_N(int nzc1,int nN1, int nzc2, int nN2);

//# cross-covariance between cluster lensing and other 2pt statistics
double cov_G_shear_cgl (double l, double delta_l,int nzs1,int nzs2, int nzc, int nN, int nzs3);
double cov_NG_shear_cgl (double l1,double l2,int nzs1,int nzs2, int nzc, int nN, int nzs3);
double cov_G_cl_cgl (double l, double delta_l,int nzl1,int nzl2, int nzc, int nN, int nzs);
double cov_NG_cl_cgl (double l1,double l2,int nzl1,int nzl2, int nzc, int nN, int nzs);
double cov_G_ggl_cgl (double l, double delta_l,int nzl,int nzs1, int nzc, int nN, int nzs2);
double cov_NG_ggl_cgl (double l1, double l2, int nzl,int nzs1, int nzc, int nN, int nzs2);

//# cross-covariance between cluster counts and 2pt statistics
double cov_shear_N (double l,int nzs1,int nzs2, int nzc, int nN);
double cov_cl_N (double l,int nzl1,int nzl2, int nzc, int nN);
double cov_ggl_N (double l,int nzl,int nzs, int nzc, int nN);

//########################
//#cluster covariance building blocks
//add integration routine for sigmaij_M


double int_I04_mmcm (double lgM, void *para){
  double *array = (double *) para;
  double u,x1,x2, c, m = exp(lgM), a= array[0];
  c = conc(m,a);
  u = pow(u_nfw_c(c,array[3],m,a),2.0)*u_nfw_c(c,array[4],m,a);
  x1 = (lgM_obs(Cluster.N_min[(int)array[1]],a)-lgM)/(sqrt(2.0)*scatter_lgM_obs(Cluster.N_min[(int)array[1]],a));
  x2 = (lgM_obs(Cluster.N_max[(int)array[1]],a)-lgM)/(sqrt(2.0)*scatter_lgM_obs(Cluster.N_max[(int)array[1]],a));
  return massfunc(m,a)*m*pow(m/(cosmology.rho_crit*cosmology.Omega_m),3.0)*u*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1));
}

double int_I04_cmcm (double lgM, void *para){
  double *array = (double *) para;
  double u,x1,x2, c, m = exp(lgM), a= array[0];
  c = conc(m,a);
  u = u_nfw_c(c,array[3],m,a)*u_nfw_c(c,array[4],m,a);
  x1 = (lgM_obs(Cluster.N_min[(int)array[1]],a)-lgM)/(sqrt(2.0)*scatter_lgM_obs(Cluster.N_min[(int)array[1]],a));
  x2 = (lgM_obs(Cluster.N_max[(int)array[1]],a)-lgM)/(sqrt(2.0)*scatter_lgM_obs(Cluster.N_max[(int)array[1]],a));
  return massfunc(m,a)*m*pow(m/(cosmology.rho_crit*cosmology.Omega_m),2.0)*u*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1));
}
double tri_1h_cmcm (double k1, double k2, double a,int nzc,int nN1, int nN2){//cmcm 1-halo term
  if (nN1 != nN2){return 0.;}// two different clusters, hence no 1-halo term
  double array[5];
  array[0]=a;
  array[1]=(double) nN1;
  array[2]=(double) nzc;
  array[3]=k1;
  array[4]=k2;
  return int_gsl_integrate_medium_precision(int_I04_cmcm, (void*)array, lgMmin(a,nN1), lgMmax(a,nN1),NULL,1000)/pow(int_gsl_integrate_medium_precision(int_n_Mobs, (void*)array, lgMmin(a,nN1),lgMmax(a,nN1),NULL,1000),2.0);
}
double tri_1h_mmcm(double k1, double k2, double a,int nzc,int nN){ //mmcm 1-halo term
  double array[5];
  array[0]=a;
  array[1]=(double) nN;
  array[2]=(double) nzc;
  array[3]=k1;
  array[4]=k2;
  return int_gsl_integrate_medium_precision(int_I04_mmcm, (void*)array, lgMmin(a,nN), lgMmax(a,nN),NULL,1000)/int_gsl_integrate_medium_precision(int_n_Mobs, (void*)array, lgMmin(a,nN),lgMmax(a,nN),NULL,1000);
//  return tri_1h_cov(k1,k2,a);
}

double int_I12_cm (double lgM, void *para){
  double x1,x2,m,a,u_m;
  double *array = (double *) para;
  m = exp(lgM); a = array[0];
  x1 = (lgM_obs(Cluster.N_min[(int)array[1]],a)-lgM)/(sqrt(2.0)*scatter_lgM_obs(Cluster.N_min[(int)array[1]],a));
  x2 = (lgM_obs(Cluster.N_max[(int)array[1]],a)-lgM)/(sqrt(2.0)*scatter_lgM_obs(Cluster.N_max[(int)array[1]],a));
  u_m = m/(cosmology.rho_crit*cosmology.Omega_m)*u_nfw_c(conc(m,a),array[2],m,a); //density profile
  return m*massfunc(m,a)*B1(m,a)*u_m*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1))*P_cm_offcentering(array[2],m,a);
}

double I12_SSC_cm(double k, double a, int nz, int nN){
  if (a < 1./1./(1+tomo.cluster_zmax[nz]) || a > 1./1./(1+tomo.cluster_zmin[nz])){return 0.;}
  double array[3] = {a,(double) nN,k};
  return int_gsl_integrate_medium_precision(int_I12_cm, (void*)array, lgMmin(a,nN), lgMmax(a,nN),NULL,1000)/int_gsl_integrate_medium_precision(int_n_Mobs, (void*)array, lgMmin(a,nN),lgMmax(a,nN),NULL,1000);
}
double delP_SSC_cm(double k, double a, double nz, double nN){
  double b= b_cluster((int)nz,(int)nN);
  return b*(68./21.+LD_term(k) -b)*p_2h(k,a)+ I12_SSC_cm(k,a,(int)nz,(int)nN);
}
//############ cluster cov routines ###########
double int_cov_N_N(double a, void *params)
{
  double *array = (double *) params;
 return  dchi_da(a)*dN200_dchi(a,(int)array[1])* dN200_dchi(a,(int)array[2])*b_cluster((int)array[0],(int)array[1])*b_cluster((int)array[0],(int)array[2])*survey_variance(a,survey.area/41253.0);
}

double cov_N_N(int nzc1,int nN1, int nzc2, int nN2){
  double res = 0;
  if (nzc1 == nzc2 && nN1 == nN2){res = N_N200(nzc1,nN1);}
  if (nzc1 == nzc2){
    double array[3] = {(double)nzc1, (double) nN1, (double) nN2};
    res += int_gsl_integrate_medium_precision(int_cov_N_N, (void*)array, 1./(1+tomo.cluster_zmax[nzc1]), 1./(1.+tomo.cluster_zmin[nzc1]),NULL,1000);
  }
  return res;
}


double project_tri_cgl_cgl(double a,void *params)
{
  double k1,k2,res = 0.,wa,weights;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
	weights = W_cluster(a,(int)ar[2],(int)ar[3]); //normalized redshift distribution of lensing cluster bin 1
  weights *= W_kappa(a,wa,ar[4]);//lens efficiency source bin 1
  weights *= W_cluster(a,(int)ar[2],(int)ar[5]); //normalized redshift distribution of lensing cluster bin 2
  weights *= W_kappa(a,wa,ar[6]);//lens efficiency source bin 2
  weights *= dchi_da(a);
  if (weights > 0){
    res = (b_cluster((int)ar[2],(int)ar[3])*b_cluster((int)ar[2],(int)ar[5])*tri_multih_cov(k1,k2,a)+tri_1h_cmcm(k1,k2,a,(int)ar[2],(int)ar[3],(int)ar[5]))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
//    res = (b_cluster((int)ar[2],(int)ar[3])*b_cluster((int)ar[2],(int)ar[5])*tri_matter_cov(k1,k2,a))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res += delP_SSC_cm(k1, a, ar[2], ar[3])* delP_SSC_cm(k2, a, ar[2], ar[5])*survey_variance(a,survey.area/41253.0)*pow(wa,-4.);
  }
  return res;
}
  

double cov_NG_cgl_cgl(double l1,double l2,int nzc1,int nN1, int nzs1, int nzc2, int nN2, int nzs2){
  if (nzc1 != nzc2){return 0.;}
  double array[7] = {l1,l2,(double)nzc1, (double) nN1, (double) nzs1, (double) nN2, (double) nzs2};  
  return  int_gsl_integrate_low_precision(project_tri_cgl_cgl,(void*) array, 1./(1+tomo.cluster_zmax[nzc1]),1./(1+tomo.cluster_zmin[nzc1]), NULL,1000);
}

double cov_G_cgl_cgl(double l, double delta_l,int nzc1,int nN1, int nzs1, int nzc2, int nN2, int nzs2){
  double C13, C14, C23, C24, N13 =0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_ccl_tomo(l,nzc1,nN1,nzc2, nN2);
  C24 = C_shear_tomo_nointerp(l,nzs1,nzs2);
  C14 = C_cgl_tomo_nointerp(l,nzc1,nN1,nzs2);
  C23 = C_cgl_tomo_nointerp(l,nzc2,nN2,nzs1);
  if (nzc1 == nzc2 && nN1 == nN2){N13= 1./n_N200(nzc1,nN1);}
  if (nzs1 == nzs2){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(nzs1)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+ C13*N24+N13*C24 + N13*N24+C14*C23)/((2.*l+1.)*delta_l*fsky);
  
}

double project_cov_cgl_N(double a,void *params)
{
  double k,res=0.,weights,wa;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k = ar[0]/wa;
  weights =W_cluster(a,(int)ar[1],(int)ar[2]); //normalized redshift distribution of lensing cluster bin
	weights *= W_kappa(a,wa,ar[3]); //lens efficiency
  weights *= dN200_dchi(a,(int)ar[4]);//dN/dV of number counts cluster bin
  weights *= dchi_da(a);
  if (weights > 0){
    res= delP_SSC_cm(k, a, ar[1], ar[2])*b_cluster((int)ar[1],(int)ar[4])*survey_variance(a,survey.area/41253.0);
  }
  return res*weights;
}


double cov_cgl_N(double l,int nzc1,int nN1, int nzs,int nzc2, int nN2){
  double res;
  if (nzc1 != nzc2){return 0.;}
  double array[5] = {l,(double)nzc1, (double) nN1, (double) nzs, (double) nN2};
  return int_gsl_integrate_medium_precision(project_cov_cgl_N,(void*) array, 1./(1+tomo.cluster_zmax[nzc1]),1./(1+tomo.cluster_zmin[nzc1]), NULL,1000);
}

//############ X-cov routines ############

double project_cov_cl_N(double a,void *params)
{
  double k,res =0.,weights,wa;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k = ar[0]/wa;
  weights =W_gal(a,(int)ar[1]); //normalized redshift distribution of galaxy sample
	weights *= W_gal(a,(int)ar[1]);
  weights *= dN200_dchi(a,(int)ar[3]);//dN/dV of number counts cluster bin
  weights *= dchi_da(a);
  if (weights > 0){
    res= (delP_SSC(k, a)-2.*bgal_a(a,ar[1])*Pdelta(k,a))*b_cluster((int)ar[2],(int)ar[3])*survey_variance(a,survey.area/41253.0);
  }
  return res*weights;
}
double project_cov_ggl_N(double a,void *params)
{
  double k,res=0.,weights,wa;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k = ar[0]/wa;
  weights =W_gal(a,(int)ar[1]); //normalized redshift distribution of galaxy sample
	weights *= W_kappa(a,wa,ar[2]);
  weights *= dN200_dchi(a,(int)ar[4]);//dN/dV of number counts cluster bin
  weights *= dchi_da(a);
  if (weights > 0){
    res= (delP_SSC(k, a)-bgal_a(a,ar[1])*Pdelta(k,a))*b_cluster((int)ar[3],(int)ar[4])*survey_variance(a,survey.area/41253.0);
  }
  return res*weights;
}
double project_cov_shear_N(double a,void *params)
{
  double k,res=0.,weights,wa;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k = ar[0]/wa;
  weights =W_kappa(a,wa,ar[1]);
	weights *= W_kappa(a,wa,ar[2]);
  weights *= dN200_dchi(a,(int)ar[4]);//dN/dV of number counts cluster bin
  weights *= dchi_da(a);
  if (weights > 0){
    res= delP_SSC(k,a)*b_cluster((int)ar[3],(int)ar[4])*survey_variance(a,survey.area/41253.0);
  }
  return res*weights;
}

double cov_cl_N(double l,int nl1,int nl2, int nzc, int nN){
  if (nl1 != nl2){return 0.;}
  if ((tomo.clustering_zmin[nl1]+tomo.clustering_zmax[nl1])/2. > tomo.cluster_zmax[nzc]){return 0.;}
  if ((tomo.clustering_zmin[nl1]+tomo.clustering_zmax[nl1])/2. < tomo.cluster_zmin[nzc]){return 0.;}
  double array[4] = {l,(double)nl1, (double) nzc, (double) nN};
  return int_gsl_integrate_medium_precision(project_cov_cl_N,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}
double cov_ggl_N(double l,int nzl,int nzs, int nzc, int nN){
  if ((tomo.clustering_zmin[nzl]+tomo.clustering_zmax[nzl])/2. > tomo.cluster_zmax[nzc]){return 0.;}
  if ((tomo.clustering_zmin[nzl]+tomo.clustering_zmax[nzl])/2. < tomo.cluster_zmin[nzc]){return 0.;}
  double array[5] = {l,(double)nzl,(double)nzs, (double) nzc, (double) nN};
  return int_gsl_integrate_medium_precision(project_cov_ggl_N,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}

double cov_shear_N(double l,int ns1,int ns2, int nzc, int nN){
  if (fmin(tomo.shear_zmin[ns1],tomo.shear_zmin[ns2]) < tomo.cluster_zmax[nzc]){return 0.;}
  double array[5] = {l,(double)ns1,(double)ns2, (double) nzc, (double) nN};
  return int_gsl_integrate_medium_precision(project_cov_shear_N,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}

double project_tri_shear_cgl(double a,void *params)
{
  double k1,k2,res=0,wa,weights;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
	weights = W_kappa(a,wa,ar[2]);
  weights *= W_kappa(a,wa,ar[3]);
  weights *= W_cluster(a,(int)ar[4],(int)ar[5]);
  weights *= W_kappa(a,wa,ar[6]);
  weights *= dchi_da(a);
  if (weights > 0){
   //res = (b_cluster((int)ar[4],(int)ar[5])*tri_multih_cov(k1,k2,a)+tri_1h_mmcm(k1,k2,a,(int)ar[4],(int)ar[5]))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res = (b_cluster((int)ar[4],(int)ar[5])*tri_matter_cov(k1,k2,a))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res += delP_SSC(k1, a)* delP_SSC_cm(k2, a, ar[4], ar[5])*survey_variance(a,survey.area/41253.0)*pow(wa,-4.);
  }
  return weights*res;
}


double cov_NG_shear_cgl (double l1,double l2,int nzs1,int nzs2, int nzc, int nN, int nzs3){
  if (fmin(tomo.shear_zmin[nzs1],tomo.shear_zmin[nzs2]) < tomo.cluster_zmax[nzc]){return 0.;}
  double array[7] = {l1,l2,(double)nzs1, (double) nzs2, (double) nzc, (double) nN, (double) nzs3};
  return int_gsl_integrate_medium_precision(project_tri_shear_cgl,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}


double cov_G_shear_cgl (double l, double delta_l,int nzs1,int nzs2, int nzc, int nN, int nzs3){
  double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_cgl_tomo_nointerp(l,nzc, nN, nzs1);C24 = C_shear_tomo_nointerp(l,nzs2,nzs3);
  C23 = C_cgl_tomo_nointerp(l,nzc, nN, nzs2);C14 = C_shear_tomo_nointerp(l,nzs1,nzs3);
  if (nzs1 == nzs3){N14= pow(survey.sigma_e,2.0)/(2.0*nsource(nzs1)*survey.n_gal_conversion_factor);}
  if (nzs2 == nzs3){N24=pow(survey.sigma_e,2.0)/(2.0*nsource(nzs2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+ C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23+N13*N24+N14*N23)/((2.*l+1.)*delta_l*fsky);
}
double project_tri_cl_cgl(double a,void *params)
{
  double k1,k2,res=0.,wa,weights;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
	weights = W_gal(a,ar[2]);
  weights *= W_gal(a,ar[2]);
  weights *= W_cluster(a,(int)ar[3],(int)ar[4]); //normalized redshift distribution of lensing cluster bin
  weights *= W_kappa(a,wa,ar[5]);
  weights *= dchi_da(a);
  if (weights > 0){
    //res = (b_cluster((int)ar[3],(int)ar[4])*tri_multih_cov(k1,k2,a)+tri_1h_mmcm(k1,k2,a,(int)ar[3],(int)ar[4]))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res = (b_cluster((int)ar[3],(int)ar[4])*tri_matter_cov(k1,k2,a))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
  res += (delP_SSC(k1, a)-2.*bgal_a(a,ar[2])*Pdelta(k1,a))* delP_SSC_cm(k2, a, ar[3], ar[4])*survey_variance(a,survey.area/41253.0)*pow(wa,-4.);
  }
  return weights*res;
}


double cov_NG_cl_cgl (double l1,double l2, int nzl1,int nzl2, int nzc, int nN,int nzs){
  if (nzl1 != nzl2){return 0.;}
  if ((tomo.clustering_zmin[nzl1]+tomo.clustering_zmax[nzl1])/2. > tomo.cluster_zmax[nzc]){return 0.;}
  if ((tomo.clustering_zmin[nzl1]+tomo.clustering_zmax[nzl1])/2. < tomo.cluster_zmin[nzc]){return 0.;}
  double array[6] = {l1,l2,(double)nzl1, (double) nzc, (double) nN, (double) nzs};
  return int_gsl_integrate_medium_precision(project_tri_cl_cgl,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}

double cov_G_cl_cgl (double l, double delta_l,int nzl1,int nzl2, int nzc, int nN,int nzs){
  double C13, C14, C23, C24;
   double fsky = survey.area/41253.0;
  C13 = C_cl_g_tomo(l,nzc, nN, nzl1);C24 = C_gl_tomo_nointerp(l,nzl2,nzs);
  C23 = C_cl_g_tomo(l,nzc, nN, nzl2);C14 = C_gl_tomo_nointerp(l,nzl1,nzs);

  return (C13*C24+  C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double cov_G_ggl_cgl (double l, double delta_l,int nzl,int nzs1, int nzc, int nN,int nzs2){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_cl_g_tomo(l,nzc,nN,nzl);C24 = C_shear_tomo_nointerp(l,nzs1,nzs2);
  C14 = C_gl_tomo_nointerp(l,nzl,nzs2);C23 = C_cgl_tomo_nointerp(l,nzc,nN,nzs1);
  if (nzs1 == nzs2){N24=pow(survey.sigma_e,2.0)/(2.0*nsource(nzs1)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+ C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23+N13*N24+N14*N23)/((2.*l+1.)*delta_l*fsky);
  
}
double project_tri_ggl_cgl(double a,void *params)
{
  double k1,k2,res=0.,wa,weights;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
	weights = W_gal(a,ar[2]);
  weights *= W_kappa(a,wa,ar[3]);
  weights *= W_cluster(a,(int)ar[4],(int)ar[5]); //normalized redshift distribution of lensing cluster bin
  weights *= W_kappa(a,wa,ar[6]);
  weights *= dchi_da(a);
  if (weights > 0){
    //res = (b_cluster((int)ar[4],(int)ar[5])*tri_multih_cov(k1,k2,a)+tri_1h_mmcm(k1,k2,a,(int)ar[4],(int)ar[5]))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res = (b_cluster((int)ar[4],(int)ar[5])*tri_matter_cov(k1,k2,a))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res += (delP_SSC(k1, a)-bgal_a(a,ar[2])*Pdelta(k1,a))*delP_SSC_cm(k2, a, ar[4], ar[5])*survey_variance(a,survey.area/41253.0)*pow(wa,-4.);
  }
  return weights*res;
}


double cov_NG_ggl_cgl (double l1,double l2,int nzl,int nzs1, int nzc, int nN,int nzs2){
  if ((tomo.clustering_zmin[nzl]+tomo.clustering_zmax[nzl])/2. > tomo.cluster_zmax[nzc]){return 0.;}
  if ((tomo.clustering_zmin[nzl]+tomo.clustering_zmax[nzl])/2. < tomo.cluster_zmin[nzc]){return 0.;}
  double array[7] = {l1,l2,(double)nzl, (double) nzs1, (double) nzc, (double) nN, (double) nzs2};
  return int_gsl_integrate_medium_precision(project_tri_ggl_cgl,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}