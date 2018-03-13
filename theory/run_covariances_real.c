void run_cov_shear_shear_real(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm1, int pm2, int start);
void run_cov_clustering_real(char *OUTFILE, char *PATH, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start);
void run_cov_ggl_real(char *OUTFILE, char *PATH, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start);
void run_cov_ggl_shear_real(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm, int start);
void run_cov_clustering_shear_real(char *OUTFILE, char *PATH, double *theta, double *dtheta,int Ntheta, int n1, int n2, int pm, int start);
void run_cov_clustering_ggl_real(char *OUTFILE, char *PATH, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start);

void run_cov_shear_shear_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int pm1, int pm2, int start);
void run_cov_ggl_real_bin(char *OUTFILE, char *PATH, double *theta,int Ntheta, int n1, int n2, int start);
void run_cov_clustering_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int start);
void run_cov_ggl_shear_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int pm, int start);
void run_cov_clustering_shear_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int pm, int start);
void run_cov_clustering_ggl_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int start);

int NG = 1;

void run_cov_shear_shear_real(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm1, int pm2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear1 = %d (%d,%d)\n", n1,z1,z2);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear2 = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  printf("%s\n",filename);
  F1 =fopen(filename,"w");
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.; c_g = 0.;
      c_g =  cov_G_shear_shear_real(theta[nl1],theta[nl2],dtheta[nl1],z1,z2,z3,z4,pm1,pm2);
      if(NG){ c_ng = cov_NG_shear_shear_real(theta[nl1],theta[nl2],z1,z2,z3,z4,pm1,pm2);}
      if(pm1==1 && pm2==1) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*n1+nl1,Ntheta*(n2)+nl2, theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm1==0 && pm2==0) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm1==1 && pm2==0) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*n1+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);  
    }
  }
  fclose(F1);
}


void run_cov_clustering_real(char *OUTFILE, char *PATH, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  //  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_1 = %d, N_cl_2 = %d\n", n1,n2);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.; 
      c_g = cov_G_cl_cl_real(theta[nl1],theta[nl2],dtheta[nl1],z1,z2,z3,z4);
      if (z1 == z3){
          if (NG){c_ng = cov_NG_cl_cl_real(theta[nl1],theta[nl2],z1,z2,z3,z4);}
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2,theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}



void run_cov_ggl_real(char *OUTFILE, char *PATH, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  zl1 = ZL(n1); zs1 = ZS(n1);
  printf("\nN_tomo_1 = %d (%d, %d)\n", n1,zl1,zs1);
  zl2 = ZL(n2); zs2 = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl2,zs2);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.;
      if(NG && zl1 == zl2 && test_zoverlap_cov(zl1,zs1)*test_zoverlap_cov(zl2,zs2)){c_ng = cov_NG_gl_gl_real(theta[nl1],theta[nl2],zl1,zs1,zl2,zs2);}
      c_g =  cov_G_gl_gl_real(theta[nl1],theta[nl2],dtheta[nl1],zl1,zs1,zl2,zs2);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);   
    }
  }
  fclose(F1);
}


void run_cov_ggl_shear_real(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm, int start)
{
  int zl,zs,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w"); 
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (test_zoverlap_cov(zl,zs)*test_zoverlap_cov(zl,z3)*test_zoverlap_cov(zl,z4) && NG){
        c_ng = cov_NG_gl_shear_real(theta[nl1],theta[nl2],zl,zs,z3,z4,pm);
      }
      c_g = cov_G_gl_shear_real(theta[nl1],theta[nl2],dtheta[nl1],zl,zs,z3,z4,pm);
      if(pm==1) fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(n2)+nl2,theta[nl1],theta[nl2],zl,zs,z3,z4,c_g,c_ng);
      if(pm==0) fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_shear_real(char *OUTFILE, char *PATH, double *theta, double *dtheta,int Ntheta, int n1, int n2, int pm, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){
        c_ng = cov_NG_cl_shear_real(theta[nl1],theta[nl2],z1,z2,z3,z4,pm);
      }
      c_g =  cov_G_cl_shear_real(theta[nl1],theta[nl2],dtheta[nl1],z1,z2,z3,z4,pm);
      if(pm==1)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(n2)+nl2,theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm==0)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2, theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
fclose(F1);
}



void run_cov_clustering_ggl_real(char *OUTFILE, char *PATH, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start)
{
  int z1,z2,zl,zs,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.; 
      c_g =  cov_G_cl_gl_real(theta[nl1],theta[nl2],dtheta[nl1],z1,z2,zl,zs);
      if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){
        c_ng = cov_NG_cl_gl_real(theta[nl1],theta[nl2],z1,z2,zl,zs);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
  }
  fclose(F1);
}


void run_cov_clustering_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  //  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_1 = %d, N_cl_2 = %d\n", n1,n2);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
  	double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
	  double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g = cov_G_cl_cl_real_rebin(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4);
      if (z1 == z3){
          if (NG){c_ng = cov_NG_cl_cl_real(t1,t2,z1,z2,z3,z4);}
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2,theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ggl_real_bin(char *OUTFILE, char *PATH, double *theta,int Ntheta, int n1, int n2, int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  zl1 = ZL(n1); zs1 = ZS(n1);
  printf("\nN_tomo_1 = %d (%d, %d)\n", n1,zl1,zs1);
  zl2 = ZL(n2); zs2 = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl2,zs2);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
  	double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
	  double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.;
      if(NG && zl1 == zl2 && test_zoverlap_cov(zl1,zs1)*test_zoverlap_cov(zl2,zs2)){c_ng = cov_NG_gl_gl_real(t1,t2,zl1,zs1,zl2,zs2);}
      c_g =  cov_G_gl_gl_real_rebin(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],zl1,zs1,zl2,zs2);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);   
    }
  }
  fclose(F1);
}
void run_cov_shear_shear_real_bin(char *OUTFILE, char *PATH,double *theta, int Ntheta, int n1, int n2, int pm1, int pm2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g,sn;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear1 = %d (%d,%d)\n", n1,z1,z2);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear2 = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
  	double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
	  double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;sn = 0.;
      if(NG){ c_ng = cov_NG_shear_shear_real(t1,t2,z1,z2,z3,z4,pm1,pm2);}
      c_g = cov_G_shear_rebin(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm1,pm2);
      if(pm1==1 && pm2==1) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*n1+nl1,Ntheta*(n2)+nl2, theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm1==0 && pm2==0) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm1==1 && pm2==0) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*n1+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);  

    }
  }
  fclose(F1);
}

void run_cov_ggl_shear_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int pm, int start)
{
  int zl,zs,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w"); 
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
  	double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
	  double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      if (test_zoverlap_cov(zl,zs)*test_zoverlap_cov(zl,z3)*test_zoverlap_cov(zl,z4) && NG){ c_ng = cov_NG_gl_shear_real(t1,t2,zl,zs,z3,z4,pm);}
      c_g = cov_G_gl_shear_real_rebin(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],zl,zs,z3,z4,pm);
      if(pm==1) fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(n2)+nl2,theta[nl1],theta[nl2],zl,zs,z3,z4,c_g,c_ng);
      if(pm==0) fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}
void run_cov_clustering_shear_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int pm, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
  	double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
	  double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real(t1,t2,z1,z2,z3,z4,pm);}
      c_g =  cov_G_cl_shear_real_rebin(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);
      if(pm==1)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(n2)+nl2,theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm==0)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2, theta[nl1],theta[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
fclose(F1);
}

void run_cov_clustering_ggl_real_bin(char *OUTFILE, char *PATH, double *theta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,zl,zs,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);
  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
  	double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
	  double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g =  cov_G_cl_gl_real_rebin(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);
      if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_real(t1,t2,z1,z2,zl,zs);}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
  }
  fclose(F1);
}

