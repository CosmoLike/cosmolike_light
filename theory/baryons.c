#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>


#include "../../theory/basics.c"
#include "../../emu13/emu.c"
#include "../../theory/parameters.c"
#include "../../theory/cosmo3D.c"
#include "../../theory/redshift.c"



double fraction_Pdelta_baryon_from_sims_DES(double k,double z)
{
  int i;
  double val,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,kintern;
  static int READ_TABLE=0;
  static double **table;
  FILE *F;
  int baryon_zbin=11;
  int baryon_kbin=100;
  double dz=(redshift.shear_zdistrpar_zmax)/(10.0);
  static double dk=(10.0-.3)/(99.);
  val = interpol2d(AGN_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0);
  return val; 
}

//   if (READ_TABLE==0){
//     if (strcmp(pdeltaparams.baryons,"AGN_DESdepth")==0) P_type = 0;
//     if (strcmp(pdeltaparams.baryons,"NOSN_DESdepth")==0) P_type = 1;
//     if (strcmp(pdeltaparams.baryons,"NOSN_NOZCOOL_DESdepth")==0) P_type = 2;
//     if (strcmp(pdeltaparams.baryons,"NOZCOOL_DESdepth")==0) P_type = 3;
//     if (strcmp(pdeltaparams.baryons,"REF_DESdepth")==0) P_type = 4;
//     if (strcmp(pdeltaparams.baryons,"WDENS_DESdepth") ==0) P_type = 5;
//     if (strcmp(pdeltaparams.baryons,"DBLIMFV1618_DESdepth")  ==0) P_type = 6;
//     if (strcmp(pdeltaparams.baryons,"WML4_DESdepth")==0) P_type = 7;
//     if (strcmp(pdeltaparams.baryons,"WML1V848_DESdepth")==0) P_type = 8;
//     if (strcmp(pdeltaparams.baryons,"AD_DESdepth")==0) P_type = 9;
//     if (strcmp(pdeltaparams.baryons,"CX_DESdepth")==0) P_type = 10;
//     if (strcmp(pdeltaparams.baryons,"CW_DESdepth")==0) P_type = 11;
//     if (strcmp(pdeltaparams.baryons,"A_DESdepth")==0) P_type = 12;
//     if (strcmp(pdeltaparams.baryons,"CSF_DESdepth")==0) P_type = 13;
//   }
//   switch (P_type){
//     case 0:  val = interpol2d(AGN_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 1:  val = interpol2d(NOSN_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 2:  val = interpol2d(NOSN_NOZCOOL_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 3:  val = interpol2d(NOZCOOL_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 4:  val = interpol2d(REF_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 5:  val = interpol2d(WDENS_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 6:  val = interpol2d(DBLIMFV1618_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 7:  val = interpol2d(WML4_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 8:  val = interpol2d(WML1V848_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 9:  val = interpol2d(AD_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 10:  val = interpol2d(CX_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 11:  val = interpol2d(CW_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 12:  val = interpol2d(A_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 13:  val = interpol2d(CSF_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;

//     default: 
//             printf("baryons.c: %s baryonic scenario not defined\n",pdeltaparams.baryons);

//             break;
//     }
//     kintern=k/cosmology.coverH0;
 
//  return val; 
// }







//    if (strcmp(pdeltaparams.baryons,"AGN_LSSTdepth")==0) P_type = 14;
//     if (strcmp(pdeltaparams.baryons,"NOSN_LSSTdepth")==0) P_type = 15;
//     if (strcmp(pdeltaparams.baryons,"NOSN_NOZCOOL_LSSTdepth")==0) P_type = 16;
//     if (strcmp(pdeltaparams.baryons,"NOZCOOL_LSSTdepth")==0) P_type = 17;
//     if (strcmp(pdeltaparams.baryons,"REF_LSSTdepth")==0) P_type = 18;
//     if (strcmp(pdeltaparams.baryons,"WDENS_LSSTdepth") ==0) P_type = 19;
//     if (strcmp(pdeltaparams.baryons,"DBLIMFV1618_LSSTdepth")  ==0) P_type = 20;
//     if (strcmp(pdeltaparams.baryons,"WML4_LSSTdepth")==0) P_type = 21;
//     if (strcmp(pdeltaparams.baryons,"WML1V848_LSSTdepth")==0) P_type = 22;
//     if (strcmp(pdeltaparams.baryons,"AD_LSSTdepth")==0) P_type = 23;
//     if (strcmp(pdeltaparams.baryons,"CX_LSSTdepth")==0) P_type = 24;
//     if (strcmp(pdeltaparams.baryons,"CW_LSSTdepth")==0) P_type = 25;
//   if (table!=0) free_double_matrix(table, 0, OWLS_kbin-1, 0, OWLS_zbin-1);
//     table   = create_double_matrix(0, OWLS_kbin-1, 0, OWLS_zbin-1);
//     printf("%s\n",file.DATA_FILE);  
//     F=fopen(file.DATA_FILE,"r");
//       for(i=0;i<OWLS_kbin;i++){
// 	fscanf(F,"%le %le %le %le %le %le %le %le %le %le %le\n",&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9,&a10,&a11);
// //	printf("%le %le %le %le %le %le %le %le %le %le %le\n",a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11);
// 	table[i][0]=a1;
// 	table[i][1]=a2;
// 	table[i][2]=a3;
// 	table[i][3]=a4;
// 	table[i][4]=a5;
// 	table[i][5]=a6;
// 	table[i][6]=a7;
// 	table[i][7]=a8;
// 	table[i][8]=a9;
// 	table[i][9]=a10;
// 	table[i][10]=a11;
//       }
//       READ_TABLE=1;
//   }
//  kintern=k/cosmology.coverH0;
//  val = interpol2d(table, OWLS_kbin, 0.3, 10.0, dk, kintern,OWLS_zbin, 0.0, 2.0, dz, z, 0.0, 0.0);
//  return val; 
// }





