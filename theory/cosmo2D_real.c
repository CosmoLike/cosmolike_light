double w_tomo_exact(int nt, int ni, int nj);// galaxy clustering tomography 2PCF galaxies in bins ni, nj, computed as sum P_l(cos(like.theta[nt])*C(l,ni,nj)
double xip_tomo_exact(int nt, int ni, int nj);//xi_+ computed as sum P_l(cos(like.theta[nt])*C(l,ni,nj)

double w_clustering_tomo(double theta, int ni, int nj); // galaxy clustering tomography 2PCF galaxies in bins ni, nj
double w_clustering_HOD(double theta, int ni); // galaxy clustering 2PCF galaxies in bin ni using HOD model

double w_gamma_t_tomo(double theta,int ni, int nj); //G-G lensing, lens bin ni, source bin nj, including IA contamination if like.IA = 3
double w_gamma_t_HOD_tomo(double theta,int ni, int nj); //G-G lensing with HOD model, lens bin ni, source bin nj

double xi_pm_tomo(int pm, double theta, int ni, int nj); //shear tomography correlation functions, including IA contamination if like.IA = 3
double xi_pm_rebin(int pm, double thetamin_i, double thetamax_i, int ni,int nj);//xi_pm averaged over large bins

/**************** these routines are only used internally ************/
typedef  double (*C_tomo_pointer)(double l, int n1, int n2);
void twopoint_via_hankel(double **xi, double *logthetamin, double *logthetamax, C_tomo_pointer C_tomo, int ni, int nj, int N_Bessel);
void xipm_via_hankel(double **xi, double *logthetamin, double *logthetamax,  C_tomo_pointer C_tomo,int ni, int nj);
// simple wrapper to make non-tomography C(l) routies fit into this pointer scheme
double C_cl_HOD_wrapper(double l,int ni, int nj){
	return C_cl_HOD(l,ni);
}
/*************** look-up tables for angular correlation functions ***************/
/******************** all angles in radian! *******************/

double w_tomo_exact(int nt, int ni, int nj){

  static int LMAX = 100000;
  static int NTHETA = 0;
  static double ** Pl =0;
  static double *Cl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  int i,l,nz;
  if (like.theta ==NULL){
    printf("cosmo2D_real.c:w_tomo_exact: like.theta not initialized\nEXIT\n");
    exit(1);
  }
  if (ni != nj){
    printf("cosmo2D_real.c:w_tomo_exact: ni != nj tomography not supported\nEXIT\n");
    exit(1);    
  }
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Cl = create_double_vector(0,LMAX-1);
    w_vec = create_double_vector(0,tomo.clustering_Nbin*like.Ntheta-1);
    NTHETA = like.Ntheta;
    char Pl_file[200];
    char Pl_file2[200];
    sprintf(Pl_file,"./aux/w_Pl_lmax%d_tmin%.1f_tmax%.1f_Nt%d",LMAX,like.vtmin/constants.arcmin,like.vtmax/constants.arcmin,NTHETA);
    sprintf(Pl_file2,"/home/teifler/Dropbox/cosmolike/top-level/des_mpp/aux/w_Pl_lmax%d_tmin%.1f_tmax%.1f_Nt%d",LMAX,like.vtmin/constants.arcmin,like.vtmax/constants.arcmin,NTHETA);
    FILE *f;
    if ((f = fopen(Pl_file, "r"))){
    	printf("reading Legendre coefficients from file %s\n",Pl_file);
	    for (i =0; i < NTHETA; i++){
    		for (int l = 0; l < LMAX; l++){
    			int j;
    			fscanf(f,"%d %d %le",&j,&j,&Pl[i][l]);
    		}
    	}
    	fclose(f);
    }
    // quick hack to make this work at JPL clusters
    //TO DO: PROPER PATH HANDLING
    else if ((f = fopen(Pl_file2, "r"))){
      printf("reading Legendre coefficients from file %s\n",Pl_file2);
      for (i =0; i < NTHETA; i++){
        for (int l = 0; l < LMAX; l++){
          int j;
          fscanf(f,"%d %d %le",&j,&j,&Pl[i][l]);
        }
      }
      fclose(f);
    }
    else{
    	for (i = 0; i<NTHETA; i ++){
      		printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
      		for (int l = 0; l < LMAX; l ++){
        		Pl[i][l] = (2.*l+1)/(4.*M_PI)*gsl_sf_legendre_Pl(l,cos(like.theta[i]));
      		}
    	}
    	printf("Legendre coefficients tabulated\nWrite to file %s\n",Pl_file);
	    if ((f = fopen(Pl_file,"w"))){
    		for (i =0; i < NTHETA; i++){
    			for (int l = 0; l < LMAX; l++){
    				fprintf(f,"%d %d %le\n",i,l,Pl[i][l]);
    			}	
		    }
		    fclose(f);
    	}
    }
  }
  if (recompute_clustering(C,G,N,ni,nj)){
    for (nz = 0; nz <tomo.clustering_Nbin; nz ++){
      for (l = 1; l < LMAX; l++){
        if (l < 20){Cl[l]=C_cl_tomo_nointerp(l,nz,nz);}
        else Cl[l]=C_cl_tomo(1.0*l,nz,nz);
      }
      for (i = 0; i < NTHETA; i++){
        w_vec[nz*like.Ntheta+i] =0;
        for (l = 1; l < LMAX; l++){
          w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
        }
      }
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[ni*like.Ntheta+nt];  
}

double xip_tomo_exact(int nt, int ni, int nj){
  static int LMAX = 20000;
  static int NTHETA = 0;
  static double ** Pl =0;
  static double *Cl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  int i,l,nz;
  if (like.theta ==NULL){
    printf("cosmo2D_real.c:xip_tomo_exact: like.theta not initialized\nEXIT\n");
    exit(1);
  }
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Cl = create_double_vector(0,LMAX-1);
    w_vec = create_double_vector(0,tomo.shear_Npowerspectra*like.Ntheta-1);
    NTHETA = like.Ntheta;
    double x,pref =1./(4.*M_PI);
    for (i = 0; i<NTHETA; i ++){
      x = cos(like.theta[i]);
      for (int l = 0; l < LMAX; l ++){
        Pl[i][l] = (2.*l+1)*pref*gsl_sf_legendre_Pl(l,x);
      }
    }
  }
  if (NTHETA != like.Ntheta){
    free_double_matrix (Pl,0, NTHETA-1, 0, LMAX-1);
    free_double_vector(w_vec,0, tomo.shear_Npowerspectra*NTHETA-1);
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    w_vec = create_double_vector(0,tomo.shear_Npowerspectra*like.Ntheta-1);
    NTHETA = like.Ntheta;
    double x,pref =1./(4.*M_PI);
    for (i = 0; i<NTHETA; i++){
      x = cos(like.theta[i]);
      for (l = 0; l < LMAX; l++){
        Pl[i][l] = (2.*l+1)*pref*gsl_sf_legendre_Pl(l,x);
      }
    }
  }
  if (recompute_shear(C,N)){
    for (nz = 0; nz <tomo.shear_Npowerspectra; nz ++){
      for (l = 0; l < LMAX; l++){
        Cl[l]=C_shear_tomo(1.0*l,Z1(nz),Z2(nz));
      //Cl[l]=C_shear_tomo_nointerp(1.0*l,Z1(nz),Z2(nz));
      }
      for (i = 0; i < NTHETA; i++){
        w_vec[nz*like.Ntheta+i] =0;
        for (l = 0; l < LMAX; l++){
          w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  return w_vec[N_shear(ni,nj)*like.Ntheta+nt];  
}
double int_for_w(double l, void *params){
  double *ar = (double *) params;
  int n1 = (int) ar[1];
  int n2 = (int) ar[2];
  return C_cl_tomo(l,n1,n2)*l*gsl_sf_bessel_J0(l*ar[0]);
}

double w_clustering_tomo(double theta, int ni, int nj) // galaxy clustering tomography 2PCF galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  if (recompute_clustering(C,G,N,ni,nj)){
    double **tab;
    int i, j,k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table   = create_double_matrix(0, tomo.clustering_Nbin*tomo.clustering_Nbin-1, 0, Ntable.N_thetaH-1);    
    for (i = 0; i < tomo.clustering_Nbin; i++){
      for (j= i; j < tomo.clustering_Nbin; j++){
        twopoint_via_hankel(tab, &logthetamin, &logthetamax,&C_cl_tomo, i,j,0);
        dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
        for (k = 0; k < Ntable.N_thetaH; k++){
          table[i*tomo.clustering_Nbin+j][k] = tab[0][k];
          table[j*tomo.clustering_Nbin+i][k] = tab[0][k];
        }
      }
    }
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);   
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return interpol(table[ni*tomo.clustering_Nbin+nj], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
}

double w_gamma_t_tomo(double theta,int ni, int nj) //G-G lensing, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res =0.;
  if (recompute_ggl(C,G,N,ni)){
  	if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: w_gamma_t_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}
  	C_tomo_pointer C_gl_pointer = &C_gl_tomo;
  	if (like.IA ==3 || like.IA ==4) C_gl_pointer = &C_ggl_IA_tab;
    double **tab;
    int i, k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table  = create_double_matrix(0, tomo.ggl_Npowerspectra, 0, Ntable.N_thetaH-1);
    for (i = 0; i <tomo.ggl_Npowerspectra; i++){
    	twopoint_via_hankel(tab, &logthetamin, &logthetamax,C_gl_pointer, ZL(i),ZS(i),2);
      for (k = 0; k < Ntable.N_thetaH; k++){table[i][k] = tab[0][k];}
    }
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH);
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if(test_zoverlap(ni,nj)) {
    res = interpol(table[N_ggl(ni,nj)], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);}
  return res;
}



double xi_pm_tomo(int pm, double theta, int ni, int nj) //shear tomography correlation functions
{
  static cosmopara C;
  static nuisancepara N; 
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  if (recompute_shear(C,N)){
  	if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: xi_pm_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}
  	C_tomo_pointer C_pointer = &C_shear_tomo;
  	if (like.IA ==3 || like.IA ==4) {C_pointer = &C_shear_shear_IA_tab;}
    update_cosmopara(&C); update_nuisance(&N);
    double **tab;
    int i,k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table   = create_double_matrix(0, 2*tomo.shear_Npowerspectra-1, 0, Ntable.N_thetaH-1);
    for (i = 0; i < tomo.shear_Npowerspectra; i++){
    	xipm_via_hankel(tab, &logthetamin, &logthetamax,C_pointer, Z1(i),Z2(i));
      for (k = 0; k < Ntable.N_thetaH; k++){
        table[2*i][k] = tab[0][k];
          table[2*i+1][k] = tab[1][k];
      }
    }
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
  }
  return interpol(table[2*N_shear(ni,nj)+(1-pm)/2], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);
}

double xi_pm_rebin(int pm, double thetamin_i, double thetamax_i, int ni,int nj){
  int Nsub_G = 10;
  int ii;
  double ti,dti;
  double xi = 0.;
  dti = (thetamax_i-thetamin_i)/(double)Nsub_G;
  for (ii = 0; ii < Nsub_G; ii++){
    ti = 2./3.*(pow(thetamin_i+(ii+1.)*dti,3.)-pow(thetamin_i+(ii+0.)*dti,3.))/(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));
    xi+=  xi_pm_tomo(pm, ti, ni, nj)*(pow(thetamin_i+(ii+1.)*dti,2.)-pow(thetamin_i+(ii+0.)*dti,2.));  
  }
  return xi/(pow(thetamax_i,2.)-pow(thetamin_i,2.));
}

double w_gamma_t_HOD_tomo(double theta,int ni, int nj) //G-G lensing with HOD model, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res = 0.;
  if (recompute_ggl(C,G,N,ni))
  {
    double **tab;
    int i, k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table  = create_double_matrix(0, tomo.ggl_Npowerspectra, 0, Ntable.N_thetaH-1);
    for (i = 0; i <tomo.ggl_Npowerspectra; i++){
    	twopoint_via_hankel(tab, &logthetamin, &logthetamax,&C_gl_HOD_tomo, ZL(i),ZS(i),2);
      for (k = 0; k < Ntable.N_thetaH; k++){table[i][k] = tab[0][k];}
    }
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if(test_zoverlap(ni,nj)) {
    res = interpol(table[N_ggl(ni,nj)], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);}
  return res;
}

double w_clustering_HOD(double theta, int ni) // HOD based galaxy clustering 2PCF galaxies in bin ni
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  if (recompute_clustering(C,G,N,ni,ni)){
    double **tab;
    int i,k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table   = create_double_matrix(0, tomo.clustering_Nbin-1, 0, Ntable.N_thetaH-1);
    
    for (i = 0; i < tomo.clustering_Nbin; i++){
        twopoint_via_hankel(tab, &logthetamin, &logthetamax,&C_cl_HOD_wrapper, i,i,0);
      dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
      for (k = 0; k < Ntable.N_thetaH; k++){
          table[i][k] = tab[0][k];
      }
    }
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return interpol(table[ni], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);
}

/****************** hankel transformation routine *******************/
void twopoint_via_hankel(double **xi, double *logthetamin, double *logthetamax, C_tomo_pointer C_tomo, int ni, int nj, int N_Bessel){
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  double loglmax, loglmin, dlnl, lnrc, arg[2];
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
  loglmax  = log(l_max);
  loglmin  = log(l_min);
  dlnl     = (loglmax-loglmin)/(1.0*Ntable.N_thetaH-1.);
  lnrc     = 0.5*(loglmax+loglmin);
  nc       = Ntable.N_thetaH/2+1;
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_tomo(l,ni,nj);
    
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
  arg[1] = N_Bessel;   /* order of Bessel function */
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

void xipm_via_hankel(double **xi, double *logthetamin, double *logthetamax,  C_tomo_pointer C_tomo,int ni, int nj)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -123.0, loglmin, dlnl,  lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i, count;
  lP   = fftw_malloc(Ntable.N_thetaH*sizeof(double));
  f_lP = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(Ntable.N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(Ntable.N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-123.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*Ntable.N_thetaH);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = Ntable.N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_tomo(l,ni,nj);
  }
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
  for (count=0; count<=1; count++) {
    arg[1] = (count==0 ? 0 : 4);   /* order of Bessel function */
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
      xi[count][Ntable.N_thetaH-i-1] = lP[i]/(t*2*constants.pi*Ntable.N_thetaH);
    }
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
