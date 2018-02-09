void init_analysis_choices();
void init_source_sample(char *sourcephotoz);
void init_lens_sample(char *lensphotoz);
double number_of_cluster_bins();

void init_analysis_choices()
{
  int i;
  // init std cosmology 
  cosmology.Omega_m   = 0.3156;
  cosmology.Omega_v   = 0.6844;
  cosmology.sigma_8   = 0.831;
  cosmology.n_spec    = 0.9645;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.0491685;
  cosmology.h0=0.6727;

  // init std survey
  survey.area   = 18000.0;
  survey.n_gal   = 26.0;
  survey.sigma_e   = 0.37;  
  survey.m_lim=27.0;
  survey.nlens=48.0;

  // binning for ell, redshift, richness bins
  like.Rmin_bias=10.0; 
  like.Ncl=25; // number of Fourier modes per tomographic power spectrum
  like.lmin= 20.0; // minimum fourier mode
  like.lmax= 20000.0; // maximum fourier mode (includes cluster WL)
  like.lmax_shear = 5000.0; // maximum fourier mode for cosmic shear

  
  // path to redshift histograms
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",SOURCE_ZFILE);
  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",LENS_ZFILE);

  // tomographic settings
  init_source_photoz("none");
  tomo.shear_Nbin=10; // tomographic bins for source sample
  tomo.shear_Npowerspectra=(int) (tomo.shear_Nbin*(tomo.shear_Nbin+1)/2); //derived quantity, number of tomographic power spectra

  set_equal_tomo_bins();
  // tomo.shear_zmin[0] = 
  // tomo.shear_zmin[1] =
  // tomo.shear_zmin[2] =
  // tomo.shear_zmin[3] =
  // tomo.shear_zmin[4] =
  // tomo.shear_zmin[5] =
  // tomo.shear_zmin[6] =
  // tomo.shear_zmin[7] =
  // tomo.shear_zmin[8] =
  // tomo.shear_zmin[9] =

  // tomo.shear_zmax[0] = tomo.shear_zmin[1];
  // tomo.shear_zmax[1] = tomo.shear_zmin[2];
  // tomo.shear_zmax[2] = tomo.shear_zmin[3];
  // tomo.shear_zmax[3] = tomo.shear_zmin[4];
  // tomo.shear_zmax[4] = tomo.shear_zmin[5];
  // tomo.shear_zmax[5] = tomo.shear_zmin[6];
  // tomo.shear_zmax[6] = tomo.shear_zmin[7];
  // tomo.shear_zmax[7] = tomo.shear_zmin[8];
  // tomo.shear_zmax[8] = tomo.shear_zmin[9];
  // tomo.shear_zmax[9] = 


  nuisance.bias_zphot_shear[0]=0.0;
  nuisance.bias_zphot_shear[1]=0.0;
  nuisance.bias_zphot_shear[2]=0.0;
  nuisance.bias_zphot_shear[3]=0.0;
  nuisance.bias_zphot_shear[4]=0.0;
  nuisance.bias_zphot_shear[5]=0.0;
  nuisance.bias_zphot_shear[6]=0.0;
  nuisance.bias_zphot_shear[7]=0.0;
  nuisance.bias_zphot_shear[8]=0.0;
  nuisance.bias_zphot_shear[9]=0.0;
  
  nuisance.sigma_zphot_shear[0]=0.05;
  nuisance.sigma_zphot_shear[1]=0.05;
  nuisance.sigma_zphot_shear[2]=0.05;
  nuisance.sigma_zphot_shear[3]=0.05;
  nuisance.sigma_zphot_shear[4]=0.05;
  nuisance.sigma_zphot_shear[5]=0.05;
  nuisance.sigma_zphot_shear[6]=0.05;
  nuisance.sigma_zphot_shear[7]=0.05;
  nuisance.sigma_zphot_shear[8]=0.05;
  nuisance.sigma_zphot_shear[9]=0.05;
  
  init_lens_photoz("none");
  tomo.clustering_Nbin=10; // tomographic bins for source sample
  tomo.clustering_Npowerspectra=tomo.clustering_Nbin;

  tomo.clustering_zmin[0] =
  tomo.clustering_zmin[1] =
  tomo.clustering_zmin[2] =
  tomo.clustering_zmin[3] =
  tomo.clustering_zmin[4] =
  tomo.clustering_zmin[5] =
  tomo.clustering_zmin[6] =
  tomo.clustering_zmin[7] =
  tomo.clustering_zmin[8] =
  tomo.clustering_zmin[9] =

  tomo.clustering_zmax[0] = tomo.clustering_zmin[1];
  tomo.clustering_zmax[1] = tomo.clustering_zmin[2];
  tomo.clustering_zmax[2] = tomo.clustering_zmin[3];
  tomo.clustering_zmax[3] = tomo.clustering_zmin[4];
  tomo.clustering_zmax[4] = tomo.clustering_zmin[5];
  tomo.clustering_zmax[5] = tomo.clustering_zmin[6];
  tomo.clustering_zmax[6] = tomo.clustering_zmin[7];
  tomo.clustering_zmax[7] = tomo.clustering_zmin[8];
  tomo.clustering_zmax[8] = tomo.clustering_zmin[9];
  tomo.clustering_zmax[9] = 

  nuisance.bias_zphot_clustering[0]=0.0;
  nuisance.bias_zphot_clustering[1]=0.0;
  nuisance.bias_zphot_clustering[2]=0.0;
  nuisance.bias_zphot_clustering[3]=0.0;
  nuisance.bias_zphot_clustering[4]=0.0;
  nuisance.bias_zphot_clustering[5]=0.0;
  nuisance.bias_zphot_clustering[6]=0.0;
  nuisance.bias_zphot_clustering[7]=0.0;
  nuisance.bias_zphot_clustering[8]=0.0;
  nuisance.bias_zphot_clustering[9]=0.0;

  nuisance.sigma_zphot_clustering[0]=0.03; 
  nuisance.sigma_zphot_clustering[1]=0.03; 
  nuisance.sigma_zphot_clustering[2]=0.03; 
  nuisance.sigma_zphot_clustering[3]=0.03; 
  nuisance.sigma_zphot_clustering[4]=0.03; 
  nuisance.sigma_zphot_clustering[5]=0.03; 
  nuisance.sigma_zphot_clustering[6]=0.03; 
  nuisance.sigma_zphot_clustering[7]=0.03; 
  nuisance.sigma_zphot_clustering[8]=0.03; 
  nuisance.sigma_zphot_clustering[9]=0.03; 


// setting cluster routines
  
  Cluster.N200_min = 10.;
  Cluster.N200_max = 220.;
  Cluster.N200_Nbin = 7;
  Cluster.l_max = like.lmax;
  Cluster.lbin=number_of_cluster_bins(); // cluster binning is a derived quantity from the above settings
  
  //N200-M relationship from Rykoff et al. 2012 (http://iopscience.iop.org/0004-637X/746/2/178/pdf/0004-637X_746_2_178.pdf) - Eq. B4 for \Delta = 200 mean version
  nuisance.cluster_Mobs_lgM0 = 1.72+log(1.e+14*0.7); //back to Msun/h instead of Msun/h70 normalization
  nuisance.cluster_Mobs_sigma = 0.25;
  nuisance.cluster_Mobs_alpha = 1.08;
  nuisance.cluster_Mobs_beta = 0.0;
  nuisance.cluster_Mobs_N_pivot = 60.;
  // completeness - adjusts later
  nuisance.cluster_completeness[0] = 0.9;
  nuisance.cluster_completeness[1] = 0.9;
  nuisance.cluster_completeness[2] = 0.9;
  nuisance.cluster_completeness[3] = 0.9;
  //no miscentering so far
  nuisance.cluster_centering_f0 = 1.0;
  nuisance.cluster_centering_alpha = 0;
  nuisance.cluster_centering_sigma = 0;
  nuisance.cluster_centering_M_pivot = 1.e+14;
  
  tomo.cluster_Nbin = 4; // number of cluster redshift bins
  tomo.cluster_zmin[0] = 0.2;
  tomo.cluster_zmax[0] = 0.4;
  tomo.cluster_zmin[1] = 0.4;
  tomo.cluster_zmax[1] = 0.6;
  tomo.cluster_zmin[2] = 0.6;
  tomo.cluster_zmax[2] = 0.8;
  tomo.cluster_zmin[3] = 0.8;
  tomo.cluster_zmax[3] = 1.0;
  tomo.cgl_Npowerspectra = 0;// number of cluster-lensing tomography combinations
  for (i = 0; i < tomo.cluster_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      tomo.cgl_Npowerspectra += test_zoverlap_c(i,j);
    }
  }

  //upper bin boundaries - note that bin boundaries need to be integers!
  int Nlist[8] = {14,20,30,45,70,120,Cluster.N200_max};
  Cluster.N_min[0] = Cluster.N200_min;
  Cluster.N_max[0] = Nlist[0];
  for (i = 1; i < Cluster.N200_Nbin; i++){
    Cluster.N_min[i] = Nlist[i-1];
    Cluster.N_max[i] = Nlist[i];
  }
 for (i = 0; i < Cluster.N200_Nbin; i++){
    printf ("Richness bin %d: %e - %e (%e Msun/h - %e Msun/h), N(z = 0.3) = %e, N(z = 0.7) = %e\n", i,Cluster.N_min[i],Cluster.N_max[i],exp(lgM_obs(Cluster.N_min[i], 0.75)),exp(lgM_obs(Cluster.N_max[i], 0.75)),N_N200(0,i),N_N200(2,i));
  }
  printf("Clusters set to LSST\n");
  printf("Clusters cgl_Npowerspectra=%d\n",tomo.cgl_Npowerspectra);
}


void set_equal_tomo_bins()
{
  int k,j;
  double frac, zi;
  
  zdistr_histo_1(0.1, NULL);
  int zbins =2000;
  double da = (redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.shear_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+zdistr_histo_1(zi, NULL);
  }
  
  tomo.shear_zmin[0] = redshift.shear_zdistrpar_zmin;
  tomo.shear_zmax[tomo.shear_Nbin-1] = redshift.shear_zdistrpar_zmax;
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.shear_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.shear_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.shear_zmax[k] = redshift.shear_zdistrpar_zmin+j*da;
    tomo.shear_zmin[k+1] = redshift.shear_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.shear_zmin[k],tomo.shear_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.shear_zmin[tomo.shear_Nbin-1],tomo.shear_zmax[tomo.shear_Nbin-1]);
  printf("redshift.shear_zdistrpar_zmin=%le max=%le\n",redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
}

void init_source_sample(char *sourcephotoz)
{
  if(strcmp(sourcephotoz,"none")==0) redshift.shear_photoz=0;
  if(strcmp(sourcephotoz,"voigt")==0) redshift.shear_photoz=1;
  if(strcmp(sourcephotoz,"voigt_out")==0) redshift.shear_photoz=2;
  if(strcmp(sourcephotoz,"gaussian")==0) redshift.shear_photoz=3;
  if(strcmp(sourcephotoz,"multihisto")==0) redshift.shear_photoz=4;
  if ((redshift.shear_photoz !=0) && (redshift.shear_photoz !=1) && (redshift.shear_photoz !=2) && (redshift.shear_photoz !=3) && (redshift.shear_photoz !=4)) 
  {
    printf("initialization error: init_source_sample: redshift.shear_photoz = %d not set properly!\nEXIT!\n",redshift.shear_photoz);
    exit(1);
  }
  printf("Source Sample Redshift Errors set to %s: redshift.shear_photoz=%d\n",sourcephotoz,redshift.shear_photoz);
}


void init_lens_sample(char *lensphotoz)
{
  if(strcmp(lensphotoz,"none")==0) redshift.clustering_photoz=0;
  if(strcmp(lensphotoz,"voigt")==0) redshift.clustering_photoz=1;
  if(strcmp(lensphotoz,"voigt_out")==0) redshift.clustering_photoz=2;
  if(strcmp(lensphotoz,"gaussian")==0) redshift.clustering_photoz=3;
  if(strcmp(lensphotoz,"multihisto")==0) redshift.clustering_photoz=4;  
    if ((redshift.clustering_photoz !=0) && (redshift.clustering_photoz !=1) && (redshift.clustering_photoz !=2) && (redshift.clustering_photoz !=3) && (redshift.clustering_photoz !=4)) 
  {
    printf("initialization error: init_lens_sample: redshift.clustering_photoz = %d not set properly!\nEXIT!\n",redshift.clustering_photoz);
    exit(1);
  }
  printf("Lens Sample Redshift Errors set to %s: redshift.clustering_photoz=%d\n",lensphotoz,redshift.clustering_photoz);
}


double number_of_cluster_bins(){
  double ell;
  int i,k=0;
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  for(i=0;i<like.Ncl;i++){
    ell=exp(log(like.lmin)+(i+0.5)*logdl);
    if (ell > like.lmax_shear){
      if (k==0) Cluster.l_min = ell;
      k=k+1;
    }
  } 
  return k;
}