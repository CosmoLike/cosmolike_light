#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

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
#include <gsl/gsl_deriv.h>

#include "basics.c"
#include "structs.c"
#include "parameters.c"
#include "../emu17/P_cb/emu.c"
#include "recompute.c"
#include "cosmo3D.c"
#include "redshift.c"
#include "halo.c"
#include "HOD.c"



int main(){
  set_cosmological_parameters_to_Planck_WP();
  limits.a_min = 1./4.;
  strcpy(pdeltaparams.runmode,"emu");
  double k;
  double a;
  for (k = 0.001; k < 1000; k*= 1.1){
	  printf("%e %e %e %e %e\n",k,p_lin(k*cosmology.coverH0,1.0),Pdelta(k*cosmology.coverH0,1.0),p_lin(k*cosmology.coverH0,.25),Pdelta(k*cosmology.coverH0,.25));
	}
  return 0;
}
