void init_DES();//set cosmology + survey parameters, specify redshift distribution files
void example_pdelta();
void example_Cl ();
void example_Cov();
void example_wtheta();
void print_zdistr();



void init_DES(){
  //survey area, etc.
	set_survey_parameters_to_DES(); //see survey.c for details
  //redshift histograms
  set_redshift_DES_conti(); //see cosmology.c for details
  //tomography bin parameters
  set_tomo_DES_conti(); //see cosmology.c for details
  //cosmological parameters
  set_cosmological_parameters_to_DES_mocks();
  set_redshift_YS();
  //	set_cosmological_parameters_to_WMAP_7years_BAO_SN(); //see cosmology.c for details
  return;
}
void example_pdelta(){
  double kmpc, k,p_l, p_nl, p_h1,p_h2,a;
  a = 1./1.45;
  printf("matter + HOD model galaxy spectra at z= 0.45: P_lin P_nl P_halo P_gg P_gm\n");
  for (kmpc = 0.001; kmpc <10.; kmpc *=1.5){
    k = kmpc*cosmology.coverH0;
    p_l = p_lin(k,a)*pow(cosmology.coverH0,3.0);
    p_nl = Pdelta (k,a)*pow(cosmology.coverH0,3.0);
    p_h1 = p_1h(k,a)*pow(cosmology.coverH0,3.0);
    p_h2 = p_2h(k,a)*pow(cosmology.coverH0,3.0);
    //    printf("%e %e %e %e %e %e\n",kmpc,p_l,p_nl,p_h2+p_h1,P_gg(k,a,0)*pow(cosmology.coverH0,3.0),P_gm(k,a,0)*pow(cosmology.coverH0,3.0) );
    printf("%e %e %e %e \n",kmpc,p_l,p_nl,p_h2+p_h1);
  }
}
void example_Cl (){
  double l;
  printf("galaxy-galaxy lensing: lens + source tomography: C_gl_tomo(l,z_l,z_s)\n");
  for (l = 20; l < 3000; l = l*1.5){ //example for calling angular power spectra
    printf ("%e %e %e %e %e %e\n",l, l*l*C_gl_tomo(l,1,0)/(2.*M_PI), l*l*C_gl_tomo(l,1,1)/(2.*M_PI),l*l*C_gl_tomo(l,1,2)/(2.*M_PI), l*l*C_gl_tomo(l,1,3)/(2.*M_PI), l*l*C_gl_tomo(l,1,4)/(2.*M_PI));
  }
  printf("\ngalaxy-galaxy lensing: lens tomography, background sources: C_gl_bg(l,z_l)\n");
  for (l = 20; l < 3000; l = l*1.5){ //example for calling angular power spectra
    printf ("%e %e %e %e %e %e\n",l, l*l*C_gl_bg(l,0)/(2.*M_PI), l*l*C_gl_bg(l,1)/(2.*M_PI),l*l*C_gl_bg(l,2)/(2.*M_PI), l*l*C_gl_bg(l,3)/(2.*M_PI), l*l*C_gl_bg(l,4)/(2.*M_PI));
  }
  printf("\ngalaxy clustering tomography (Limber approx): C_cl_tomo (l,z1,z2)\n");
  for (l = 20; l < 3000; l = l*1.5){ //example for calling angular power spectra
    printf ("%e %e %e %e %e %e\n",l, l*l*C_cl_tomo(l,0,0)/(2.*M_PI), l*l*C_cl_tomo(l,0,1)/(2.*M_PI),l*l*C_cl_tomo(l,0,2)/(2.*M_PI), l*l*C_cl_tomo(l,0,3)/(2.*M_PI), l*l*C_cl_tomo(l,4,4)/(2.*M_PI));
  }
}
void example_Cov(){
  printf("Examples for Gaussian Covariances\n");
  printf ("%e %e  %e %e\n",cov_G_Glensing(1.*constants.arcmin,2.*constants.arcmin,.25*constants.arcmin,0,2,0,4),cov_G_Glensing(1.*constants.arcmin,1.*constants.arcmin,.25*constants.arcmin,0,2,0,4),cov_G_Glensing(1.*constants.arcmin,2.*constants.arcmin,.25*constants.arcmin,0,4,0,4),cov_G_Glensing(1.*constants.arcmin,1.*constants.arcmin,.25*constants.arcmin,0,4,0,4));
  printf ("%e %e  %e %e\n",cov_G_clustering(1.*constants.arcmin,2.*constants.arcmin,1.*constants.arcmin,2,2,2,2),cov_G_clustering(1.*constants.arcmin,1.*constants.arcmin,1.*constants.arcmin,2,2,2,2),cov_G_clustering(100.*constants.arcmin,101.*constants.arcmin,1.*constants.arcmin,2,2,2,2),cov_G_clustering(100.*constants.arcmin,100.*constants.arcmin,1.*constants.arcmin,2,2,2,2));
  printf("Example for Non-Gaussian Covariances - this may take a while...\n");
  printf ("%e %e\n",cov_NG_Glensing(1.*constants.arcmin,2.*constants.arcmin,0,2,0,4),cov_NG_Glensing(1.*constants.arcmin,1.*constants.arcmin,0,2,0,4));
}
void example_wtheta(){ //angular 2PT functions using linear bias + non-linear matter power spectrum
  double t,l;
  printf("angular 2PT functions: theta [arcmin], clustering (zi,zj), g-g lensing tomography (zl,zs), g-g lensing (zl)");
  for (t = .6/1.2589; t <480; t=t*1.2589){
    printf ("%e %e \n",t/60.,w_clustering_tomo(t*constants.arcmin,0,0));//,w_gamma_t_tomo(t*constants.arcmin,0,1),w_gamma_t_bg(t*constants.arcmin,0));
  }
  printf("\ngalaxy clustering tomography (Limber approx): C_cl_tomo (l,z1,z2)\n");
  for (l = 2; l < 2.e+4; l = l*1.25){ //example for calling angular power spectra
    printf ("%e %e \n",l, l*(l+1)*C_cl_tomo(l,0,0)/(2.*M_PI));
  }
  
}
void print_zdistr(){
  double z;
  printf("\nredshift distributions: z n_lens(z,0) n_lens(z,1) n_source(z,2) n_source(z,3) n_source_bg(z,0)\n");
  for (z = 0.3; z< 0.6; z +=0.01){ //print redshift distributions
 		printf ("%e  %e %e %e %e %e \n",z, pf_photoz(z,0),pf_photoz(z,1), zdistr_photoz(z,2),zdistr_photoz(z,3),zdistr_photoz_bg(z,0));
	}
}
void example_w_HOD(){ //examples for angular correlation functions using HOD model
  /* HOD functions: angular correlation functions */
  printf(" galaxy clustering: w_clustering_HOD(theta,z), G-G lensing tomography: w_gamma_t_HOD_tomo(theta,zl,zs), G-G lensing, background sources: w_gamma_t_HOD_bg(theta, zl)\n");
  printf("%e %e %e\n",w_clustering_HOD(0.1*constants.arcmin,0),w_gamma_t_HOD_tomo(0.1*constants.arcmin,0,1),w_gamma_t_HOD_bg(0.1*constants.arcmin,0));
}

int main(void){
  int i;
  gsl_set_error_handler_off ();
  
  init_DES(); //set stardard cosmology + survey parameters + redshift distributions
  //  set_HOD(-1); //specify HOD paramers in this routine
  //  set_HOD(0);
  //  set_HOD(3);
  //use photo-z distribution? 0 -> off, 1-> on, input distribution binned in true redshift
  
  
  //initialize linear bias: b_1 (<z>) for each bin, code assumes redshift evolution as specified in routine bgal_z(z,i) within each bin
  for (i = 0; i < 5; i++){
    gbias.b[i][0] = 1.0;
    gbias.b[i][1] = (tomo.clustering_zmin[i] + tomo.clustering_zmax[i])/2.;
  }
  printf("Omega_m %f Omega_b %f\n", cosmology.Omega_m, cosmology.omb);
  //double k;
  //  for (k = 12.+0.5/15.; k < 16; k+=1./15.){
  //    printf("%e  %e  %e \n", pow(10.,k) ,massfunc(pow(10.,k),1./1.45)*pow(10.,2.*k)/cosmology.rho_m, B1(pow(10.,k),1./1.45));
  //  }
  // for (k = 1e-3; k<1.e+5; k*=1.1){
  //    double p = Pdelta(k, 0.99);
  //   printf("%le  %le\n",k,Delta_L(k));
  //  }
  
  /****** example routines for angular correlation functions, angular power spectra, covariances *****/
  /******    note that all angles need to be in radian! *******/
  
  //  example_wtheta(); //angular 2PT functions using linear bias + non-linear matter power spectrum
  //  example_w_HOD();  //angular 2PT functions using HOD model
  // example_Cl();     //angular power spectra
  // example_Cov();    //Covariance routines
  example_pdelta(); //linear + non-linear matter power spectrum, galaxy power spectrum using HOD
  //  print_zdistr();   //redshift distributions
	return 0;
}