#ifndef SCATTER_H
#define SCATTER_H

/* Module containing all variables and functions related to scattering simulations */

#include "jurassic.h"
#include "forwardmodel.h"
#include "misc.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Compute orthonormal basis. */
void bascoord(double *dz,
	      double *dy,
	      double *ex,
	      double *ey,
	      double *ez);

/* Compute Mie parameters. */
void bhmie(double x,
	   double n_real,
	   double n_imag,
	   double *phase,
	   double *qext,
	   double *qsca);

/* Gauss Hermite abcissas and weights. */
void gauher(double *x,
	    double *w);

/* Copy and initialize aerosol variables. */
void copy_aero(ctl_t *ctl,
	       aero_t *aero_dest,
	       aero_t *aero_src,
	       int init);

/* Get aerosol/cloud optical properties (1D). */
void get_opt_prop(ctl_t *ctl,
		  aero_t *aero);

/* Calculate optical properties with Mie theory for a log-normal mode. */
void opt_prop_mie_log(ctl_t *ctl,
		    aero_t *aero,
		    int count,
		    double *beta_ext,
		    double *beta_sca,
		    double phase[NDMAX][NTHETA]);

/* Get optical properties from external database. - New */
void opt_prop_external(ctl_t *ctl,
		      aero_t *aero,
		      int count,
		      double *beta_ext,
		      double *beta_sca,
		      double phase[NDMAX][NTHETA]);

/* Read aerosol/cloud data. */
void read_aero(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       aero_t *aero);

/* Compute scattering source. */
void srcfunc_sca(ctl_t *ctl,
		 atm_t *atm,
		 aero_t *aero,
		 double sec,
		 double *x,
		 double *dx,
		 int il,
		 double *src_sca,
		 int scattering,
     queue_t *q);

/* Compute scattering source (thermal emissions). */
void srcfunc_sca_1d(ctl_t *ctl,
		    atm_t *atm,
		    aero_t *aero,
        double sec,
		    double *x,
		    double *dx,
		    int il,
		    double *src_sca,
		    int scattering,
        queue_t *q);

/* Compute scattering source (thermal emissions). */
void srcfunc_sca_3d(ctl_t *ctl,
		    atm_t *atm,
		    aero_t *aero,
        double sec,
		    double *x,
		    double *dx,
		    int il,
		    double *src_sca,
		    int scattering,
        queue_t *q);

/* Add solar radiation to scattering source. */
void srcfunc_sca_sun(ctl_t *ctl,
		     atm_t *atm,
		     aero_t *aero,
		     double sec,
		     double *x,
		     double *dx,
		     int il,
		     double *src_sun,
         queue_t *q);

/* Compute Sun's angular coordinates. */
void suncoord(double sec,
	      double lon,
	      double lat,
	      double *azi,
	      double *sza);

/* Write particle data. */
void write_aero(const char *dirname,
		const char *filename,
		aero_t *aero);

#endif
