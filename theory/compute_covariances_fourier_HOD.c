/*#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <string.h>*/

/*#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>

#include "basics.c"
#include "structs.c"
#include "parameters.c"
#include "../emu13/emu.c"
#include "recompute.c"
#include "cosmo3D.c"
#include "redshift.c"
#include "halo.c"
#include "HOD.c"
#include "redmagic.c"
#include "cosmo2D_fourier.c"
#include "IA.c"
#include "cluster.c"
#include "BAO.c"
#include "external_prior.c"
#include "covariances_3D.c"
#include "covariances_fourier.c"
#include "covariances_cluster.c"*/
#include "compute_covariances_fourier.c"
#include "redmagic.c"
#include "covariances_fourier_HOD.c"

void run_cov_ggl_N_HOD (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_ggl_cgl_HOD (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_cl_N_HOD (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_cl_cgl_HOD (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);

void run_cov_ggl_shear_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_shear_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_ggl_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ggl_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);

void run_cov_ggl_N_HOD (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight){
        cov =cov_ggl_N(ell[nl1],zl,zs, nzc2, nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., zl, zs, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_ggl_cgl_HOD (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2, zl, zs, nzs1, nl2, nzc2, nzs3,i,j;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
    for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
      i =   like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
      j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;
        
      c_g = 0; c_ng = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight){
        c_ng = cov_NG_ggl_cgl(ell[nl1],ell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);
        if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){c_g =cov_G_ggl_cgl_HOD(ell[nl1],dell[nl1],zl,zs, nzc2, nN2,nzs3);}
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], zl, zs, nzc2, nzs3,c_g, c_ng);
    }
  }
  fclose(F1);
}

void run_cov_cl_N_HOD (char *OUTFILE, char *PATH, double *ell, double *dell,int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],N1);
      if (weight){
        cov =cov_cl_N(ell[nl1],N1,N1, nzc2, nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., N1, N1, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_cl_cgl_HOD (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2,nzc2, nzs3,i,j,nl2;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;
        
        c_g = 0; c_ng = 0.;
        weight = test_kmax(ell[nl1],N1);
        if (weight){
          c_ng = cov_NG_cl_cgl(ell[nl1],ell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){c_g =cov_G_cl_cgl_HOD(ell[nl1],dell[nl1],N1,N1, nzc2, nN2,nzs3);}
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
        //printf("%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
      }
    }
  fclose(F1);
}

void run_cov_ggl_shear_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl,zs,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w"); 
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(zl,z3)*test_zoverlap(zl,z4)){c_ng = cov_NG_gl_shear_tomo(ell[nl1],ell[nl2],zl,zs,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_gl_shear_tomo_HOD(ell[nl1],dell[nl1],zl,zs,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_shear_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],z1);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(z1,z3)*test_zoverlap(z1,z4)){
        c_ng = cov_NG_cl_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
      }
      if (nl1 == nl2){
        c_g =  cov_G_cl_shear_tomo_HOD(ell[nl1],dell[nl1],z1,z2,z3,z4);
      }
    }
    fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_ggl_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,zl,zs,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (z1 == zl){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],zl);
        if (weight){
          c_ng = cov_NG_cl_gl_tomo(ell[nl1],ell[nl2],z1,z2,zl,zs);
          if (nl1 == nl2){
            c_g =  cov_G_cl_gl_tomo_HOD(ell[nl1],dell[nl1],z1,z2,zl,zs);
          }
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_2 = %d\n", n2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (z1 == z3){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],z3);
        if (weight) {
          c_ng = cov_NG_cl_cl_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_cl_cl_tomo_HOD(ell[nl1],dell[nl1],z1,z2,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ggl_HOD(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2, weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  zl1 = ZL(n1); zs1 = ZS(n1);
  printf("\nN_tomo_1 = %d (%d, %d)\n", n1,zl1,zs1);
  zl2 = ZL(n2); zs2 = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl2,zs2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl1)*test_kmax(ell[nl2],zl2);
      if (weight && zl1 == zl2) {
        c_ng = cov_NG_gl_gl_tomo(ell[nl1],ell[nl2],zl1,zs1,zl2,zs2);
      }
      if (nl1 == nl2){
        c_g =  cov_G_gl_gl_tomo_HOD(ell[nl1],dell[nl1],zl1,zs1,zl2,zs2);
      }
      if (weight ==0 && n2 != n1){
        c_g = 0;
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);   
    }
  }
  fclose(F1);
}


int main(int argc, char** argv)
{
  int i,l,m,n,o,s,p,nl1;
  int hit=atoi(argv[1]);
  char OUTFILE[400],PATH[400];
  
  //RUN MODE setup

  init_cosmo();
  init_binning(0.1);
  init_survey("LSST");
  init_galaxies("../zdistris/zdistribution_LSST","../zdistris/zdistribution_const_comoving", "gaussian", "gaussian", "redmagic");
  init_HOD_rm();
  init_clusters();
  init_IA("none","GAMA");
//  init_priors("none","none","PhotoBAO");
  sprintf(PATH,"../top-level/covparallel/HOD_");
  int k=1;
  
  //set l-bins for shear, ggl, clustering
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  double *ell, *dell;
  ell=create_double_vector(0,like.Ncl-1);
  dell=create_double_vector(0,like.Ncl-1);
  for(i=0;i<like.Ncl;i++){
    ell[i]=exp(log(like.lmin)+(i+0.5)*logdl);
    dell[i]=exp(log(like.lmin)+(i+1)*logdl) - exp(log(like.lmin)+(i*logdl));
  } 
  //set l-bins for cluster lensing
  double *ell_Cluster, *dell_Cluster;
  logdl=(log(Cluster.l_max)-log(Cluster.l_min))/Cluster.lbin;
  ell_Cluster=create_double_vector(0,Cluster.lbin-1);
  dell_Cluster=create_double_vector(0,Cluster.lbin-1);
  for(i=0;i<Cluster.lbin;i++){
    ell_Cluster[i]=exp(log(Cluster.l_min)+(i+0.5)*logdl);
    dell_Cluster[i]=exp(log(Cluster.l_min)+(i+1)*logdl) - exp(log(Cluster.l_min)+(i*logdl));
  }

  sprintf(OUTFILE,"%s_ssss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.shear_Npowerspectra; l++){
    for (m=l;m<tomo.shear_Npowerspectra; m++){
      if(k==hit) run_cov_shear_shear(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
      //printf("%d\n",k);
    }
  }

  printf("%d\n",k);
  sprintf(OUTFILE,"%s_lsls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
    for (m=l;m<tomo.ggl_Npowerspectra; m++){
      if(k==hit) run_cov_ggl_HOD(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
    }
  }
  sprintf(OUTFILE,"%s_llll_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){ //auto bins only for now!
    for (m=l;m<tomo.clustering_Npowerspectra; m++){
      if(k==hit) run_cov_clustering_HOD(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
  //    printf("%d %d %d\n",l,m,k);
    }
  }
  sprintf(OUTFILE,"%s_llss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      if(k==hit) run_cov_clustering_shear_HOD(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
    //  printf("%d\n",k);
    }
  }
  sprintf(OUTFILE,"%s_llls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
    for (m=0;m<tomo.ggl_Npowerspectra; m++){
      if(k==hit) run_cov_clustering_ggl_HOD(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
      //printf("%d\n",k);
    }
  }
  sprintf(OUTFILE,"%s_lsss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      if(k==hit) run_cov_ggl_shear_HOD(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
    }
  }
printf("%d\n",k);

/********** cluster covariance ************/
  sprintf(OUTFILE,"%s_nn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.cluster_Nbin; l++){
    for (m=0;m<tomo.cluster_Nbin; m++){
      if(k==hit) run_cov_N_N (OUTFILE,PATH,l,m,k);
      k=k+1;
      //printf("%d\n",k);
    }
  }
  sprintf(OUTFILE,"%s_cscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.cgl_Npowerspectra; l++){
    for (m=0;m<tomo.cgl_Npowerspectra; m++){
      if(k==hit) run_cov_cgl_cgl (OUTFILE,PATH,ell_Cluster,dell_Cluster,l,m,k);
      k=k+1;
    }
  }
  sprintf(OUTFILE,"%s_csn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.cgl_Npowerspectra; l++){
    for (m=0;m<tomo.cluster_Nbin; m++){
      if(k==hit) run_cov_cgl_N (OUTFILE,PATH,ell_Cluster,dell_Cluster,l,m,k);
      k=k+1;
    }
  }
// shear X cluster
  sprintf(OUTFILE,"%s_ssn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.shear_Npowerspectra; l++){
    for (m=0;m<tomo.cluster_Nbin; m++){
      if(k==hit) run_cov_shear_N (OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
    }
  }  
  sprintf(OUTFILE,"%s_sscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.shear_Npowerspectra; l++){
    for (m=0;m<tomo.cgl_Npowerspectra; m++){
      for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
        if(k==hit) run_cov_shear_cgl (OUTFILE,PATH,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
        k=k+1;
      }
    }
  }
  // ggl X cluster
  printf("%d\n",k);
  sprintf(OUTFILE,"%s_lsn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
    for (m=0;m<tomo.cluster_Nbin; m++){
      if(k==hit) run_cov_ggl_N_HOD (OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
    }
  }
  sprintf(OUTFILE,"%s_lscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
    for (m=0;m<tomo.cgl_Npowerspectra; m++){
      for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
        if(k==hit) run_cov_ggl_cgl_HOD (OUTFILE,PATH,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
        k=k+1;
        //printf("%d\n",k);
      }
    }
  }
  // clustering X cluster
  sprintf(OUTFILE,"%s_lln_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
    for (m=0;m<tomo.cluster_Nbin; m++){
      if(k==hit) run_cov_cl_N_HOD (OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
      //printf("%d\n",k);
    }
  }
  sprintf(OUTFILE,"%s_llcs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
    for (m=0;m<tomo.cgl_Npowerspectra; m++){
      for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
        if(k==hit) run_cov_cl_cgl_HOD (OUTFILE,PATH,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
        k=k+1;
        //printf("%d\n",k);
      }
    }
  }
  
  printf("number of cov blocks for parallelization: %d\n",k-1); 
  
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}

