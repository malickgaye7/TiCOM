/********************************************************
 * conduction -- Routines to compute the conductive	    *
 *	temperature structure for a spherically 	          *
 *	symmetric planet or moon composed of an 	          *
 *	arbitrary number of layers.			                    *
 *							                                        *
 * Purpose:						                                  *
 *	Originally designed to see if the heat flux	        *
 *	generated by tidal dissipation in an ice shell	    *
 *	is consistent with a conductive temperature	        *
 *	profile.  Because the tidal heating is a 	          *
 *	function of the temperature-dependent ice 	        *
 *	viscosity, and the shell thickness, the heating	    *
 * 	should be solved at every iteration as well.	      *
 *	This code should be coupled with tidal_heating	    *
 *	code in order to do this.			                      *
 *							                                        *
 * Notes:						                                    *
 *	Based on a freestanding code, conduction.c	        *
 *							                                        *
 * Author:  James Roberts				                        *
 *							                                        *
 * Date:  13 April 2007					                        *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "tirade.h"

/* Make input file params global, so we don't have to pass
 * them around all day. */

extern long double **heating_axi_ice;
extern long double **T;

extern struct grid_info grids;
extern struct globals planet;
extern struct boundaries bc;
extern struct matprop *layers;

void conduction(char directory[], int tc_iter)
{

    long double misfit,mf_prev;	/* Misfit between calc & obs. T at surface */
    long double r_core = 1.6e5;	/* Core radius */

    int i,j,k;		/* Counter */
    int iter;		/* Iteration of T solver */
    char filename[80];			/* Test File name */

    FILE *out_file1;


    void initialize_temperature(long double **, char directory[]);
    void compute_temp(long double **, long double **, long double *, 
	    char directory[], int, int);

    /* Verify heating data from tidal module.  Can remove for final vesion. */
/*    fprintf(stderr,"r_inner %Lf r_outer %Lf \n",grids.r_inner,grids.r_outer);
    fprintf(stderr,"dth %.1Lf dr %.2Le \n",grids.dth*180.0/M_PI,grids.dr);

    sprintf(filename,"%s/iceh.dat",directory);
    out_file1 = fopen(filename,"w");
    for(i=1;i<=grids.num_theta;i++)
	for(k=1;k<=(grids.num_r-3);k++)
	    (void)fprintf(out_file1,"%.1Lf %.2Le %Le %d\n",
		grids.etheta[i]*180.0/M_PI,grids.ice_radius[k],
		heating_axi_ice[i][k],k);
    fclose(out_file1);
*/
    /* Compute temperature profile */
    if (tc_iter == 0)
	initialize_temperature(T, directory);   /* Make Initial Guess */
    iter = 0;
    compute_temp(T, heating_axi_ice, &misfit, directory, tc_iter, iter);
    //fprintf(stderr,"misfit %d = %Le, ri = %Lf\n",iter,misfit,grids.r_inner);

    /* Adjust base of ice shell and re-compute T until misfit becomes small */
    while (fabs(misfit) > 0.01  && iter < 100) {

	if (fabs(misfit) < 1.0)
	    grids.r_inner += -0.1*misfit * (grids.r_outer 
						- grids.r_inner);
	else grids.r_inner += -misfit/fabs(misfit) * 1000.0;
	if (grids.r_inner >= 0.999 * grids.r_outer) {
	    (void)fprintf(stderr,"Shell melted:  No solution\n");
	    break;
	}
	else if (grids.r_inner <= r_core) {
	    (void)fprintf(stderr,"Ocean froze:  No solution\n");
	    break;
	}
	
	/* Adjust radial grid information */
	/* r_inner has been adjusted already */

	grids.dr = (grids.r_outer - grids.r_inner) 
		    	    / ((long double)(grids.num_r-3));

	for(k=1;k<=(grids.num_r-1);k++) {
	    grids.ice_radius[k] = grids.r_outer 
		- (grids.r_outer - grids.ice_radius[k]) 
		    / (grids.r_outer - grids.ice_radius[0]) 
		    * (grids.r_outer - grids.r_inner);
		    
	    grids.eradius[k+2] = grids.ice_radius[k];

	    if (k >= 1)
		grids.radius[k+1] = 0.5*(grids.eradius[k+1] 
					    + grids.eradius[k+2]);
	}
	/* Set ice_radius[0], radius[2] to r_inner, radius[nr-1] to r_outer */
	grids.ice_radius[0] = grids.r_inner;
	grids.radius[2] = grids.r_inner;
	grids.radius[grids.num_r-1] = grids.r_outer;

	/* Need to fix eradius[2] (center of ocean) */
	grids.eradius[2] = 0.5*(grids.radius[1] + grids.radius[2]);

	mf_prev = misfit;
	iter++;
	
	compute_temp(T, heating_axi_ice, &misfit, directory, tc_iter, iter);
//	fprintf(stderr,"misfit %d = %Lf, ri = %Lf\n",iter,misfit,grids.r_inner);

/*	if (fabs(misfit) >= fabs(mf_prev)) {
	    (void)fprintf(stderr,"No convergence:  No solution\n");
	    break;
	}
*/	
    }

    /* Print new coordinates */
/*  sprintf(filename,"%s/coords1.dat",directory);
    out_file1 = fopen(filename,"w");
    for(i=0;i<=grids.num_theta-1;i++) 
	for(j=0;j<=grids.num_phi-1;j++)
	    for(k=0;k<=grids.num_r-1;k++) 
		(void)fprintf(out_file1,"%Lf %Lf %Lf\n",
			grids.theta[i],grids.phi[j],grids.radius[k]);
    fclose(out_file1);

    sprintf(filename,"%s/coords_elem1.dat",directory);
    out_file1 = fopen(filename,"w");
    for(i=0;i<=grids.num_theta;i++)
	for(j=0;j<=grids.num_phi;j++)
	    for(k=0;k<=grids.num_r;k++)
		(void)fprintf(out_file1,"%Lf %Lf %Lf\n",
			grids.etheta[i],grids.ephi[j],grids.eradius[k]);
    fclose(out_file1);

    sprintf(filename,"%s/coords_ice_elem1.dat",directory);
    out_file1 = fopen(filename,"w");
    for(i=0;i<=grids.num_theta;i++)
	for(j=0;j<=grids.num_phi;j++)
	    for(k=0;k<=(grids.num_r-2);k++)
		(void)fprintf(out_file1,"%Lf %Lf %Lf\n",
			grids.etheta[i],grids.ephi[j],grids.ice_radius[k]);
    fclose(out_file1);
*/
}


/************************************************************
 * initialize_temperature -- Makes an initial guess at the  *
 *		solution for the temperature and outputs it.          *
 *							                                            *
 * Parameters						                                    * 	
 *	T	    -- Temperature			                              *
 *	directory   -- Data directory			                      *
 ************************************************************/
void initialize_temperature(long double **T, char directory[]) {

    FILE *out_file;	    /* Output file */
    char output_file[80];   /* Full path of temperature file */
    int i,k;		    /* Counters */
    int nr,nth;	    
    long double dth,dr;	    /* deltas */
    long double thi,rk;	    /* current coords */
    long double latdep;	/* variable surface temperature */

    nr = grids.num_r - 2;   /* Ignore silicate and ocean layers */
    nth = grids.num_theta;

    dth = grids.dth;
    dr = grids.dr;

    /* Apply BCs */
    for (i=0;i<=nth;i++) {
	thi = grids.etheta[i];
	T[i][0] = bc.bottemp;
	T[i][nr] = bc.surftemp;

	/* Use Variable Surface Temperature */
	/* Ojakangas & Stevenson (1989) */

	if ( (thi >= planet.ob) && (thi <= (M_PI-planet.ob)) )
           latdep = pow(sin(thi),0.25);
        else if (thi < planet.ob)
           latdep = pow( ((thi*thi + planet.ob*planet.ob)/2.0), 0.125);
        else
           latdep = pow( (((M_PI-thi)*(M_PI-thi) + planet.ob*planet.ob)/2.0),0.125);

	T[i][nr] *= latdep;
//	(void)fprintf(stderr,"latdep %d %Lf %Lf %Lf\n",i,thi,latdep,T[i][nr]);

    }

    /* Initial Guess at T, linear in r, const in th */
    for (k=1;k<=(nr-1);k++) {
	rk = grids.ice_radius[k];
	for (i=0;i<=nth;i++) {
	    thi = grids.etheta[i];
	    T[i][k] = bc.surftemp + (bc.bottemp - bc.surftemp) 
		* (grids.r_outer - rk) 
		/ (grids.r_outer - grids.r_inner);
	}
    }

    sprintf(output_file,"%s/Tcond.initial",directory);
    out_file = fopen(output_file,"w");

    for (i=0;i<=nth;i++) {
	thi = grids.etheta[i] * 180.0/M_PI;
	for (k=0;k<=nr;k++) {
	    rk = grids.ice_radius[k]/grids.r_outer;
	    fprintf(out_file,"%.2Lf %.3Le %.2Lf\n",thi,rk,T[i][k]);
	}
    }

    fclose(out_file);

}

/************************************************************
 * compute_temp -- solve for the temperature using a 	      *
 *		finite-difference CS scheme and calculate             *
 *		the misfit at the surface.		                        *
 *							                                            *
 * Parameters						                                    *
 *	H	    -- Volumetric heating rate		                    *
 *	T	    -- Temperature			                              *
 *	misfit_ptr  -- misfit between calc and obs. Tsurf       *
 *	directory   -- data directory			                      *
 *	iter	    -- Shell thickness iteration	                *
 ************************************************************/

void compute_temp(long double **T, long double **H, long double *misfit_ptr, 
	char directory[], int tc_iter, int iter) {

    FILE *out_file;	    /* Output file */
    char output_file[80];   /* Full path of heating_file */
    int i,k;		    /* Counters */
    int gs_iter;	    /* Gauss-Seidel counter */
    int nr,nth;	    
    long double dth,dr;	    /* deltas */
    long double thi, rk;	    /* current coords */
    long double max_eps;	    /* Error in current GS iteration */
    long double flux;		    /* Heat flux at bottom */
    long double scaledflux;	    /* Core flux scaled to r at base */

    long double gauss_seidel(long double **, long double **);

    nr = grids.num_r - 2;  /* Ignore silicate and ocean layers */
    nth = grids.num_theta;

    dth = grids.dth;
    dr = grids.dr;

    gs_iter = 0;
    max_eps = 1.0;

    /* Gauss-Seidel Iteration */
    //while (gs_iter < (nr*nth) && max_eps > 0.001) {
    while (gs_iter < 1000 && max_eps > 0.001) {
	gs_iter++;
	max_eps = gauss_seidel(T,H);
    
	//(void)fprintf(stdout,"G-S Iter = %d, error = %Le\n",gs_iter,max_eps);
    }  
    //(void)fprintf(stderr,"G-S Iter = %d, error = %Le, ",gs_iter,max_eps);
 
    sprintf(output_file,"%s/Tcond.%d.%d",directory,tc_iter,iter);
    out_file = fopen(output_file,"w");

    /* Output Temp & Calculate Heat Flux at Base */

    flux = 0.0;

    for (i=0;i<=nth;i++) {
	thi = grids.etheta[i];
	for (k=0;k<=nr;k++) {
	    rk = grids.ice_radius[k]/grids.r_outer;
	    fprintf(out_file,"%Lf %Lf %Lf %Le\n",thi*180.0/M_PI,rk,T[i][k],H[i][k]);
	}

	flux += (T[i][1] - T[i][0]) * sin(thi);
/*	if(i==1 || i==19)
	    (void)fprintf(stderr,"flux %d %Lf %Lf %Lf %Lf\n",
		i,thi,T[i][1],T[i][2],flux);*/
    }

    fclose(out_file);
    flux *= 0.5 * dth;	/* 0.5 needed, Integral of sin(th) from 0 to Pi = 2 */
    (void)fprintf(stderr,"dT %Lf dr %Lf ",flux,0.5*dr);
    flux *= -planet.k / (0.5*dr);  /* only 1/2 an element btw 1st el 
					  and bottom */
    (void)fprintf(stderr,"flux %Lf\n",flux);

    scaledflux = bc.flux 
	* (layers[1].rad/grids.r_inner) * (layers[1].rad/grids.r_inner);  
    (void)fprintf(stderr,"flux bc %Lf\n",scaledflux);
    *misfit_ptr = ((flux - scaledflux) / scaledflux);
  
}


/************************************************************
 * gauss_seidel -- One iteration of the Gauss-Seidel	      *
 *		solution for the temperature.  Calculate              *
 *		the error from the previous iteration.	              *
 *							                                            *
 * Parameters						                                    *	
 *	T	    -- Temperature			                              *
 *	H	    -- Volumetric heating rate		                    *
 *							                                            *
 * Returns						                                      *
 *	max_eps	-- Error from previous iteration.	              *
 ************************************************************/
long double gauss_seidel(long double **T, long double **H) {

    long double thi, rk;			/* Current coordinates */
    long double dth, dth2, dr, delta, delta2;	/* Spatial Intervals */
    long double Tprev;			/* Temperature at prev. iter. */
    long double eps, max_eps;			/* Error */
    int i,k;				/* counters */
    int nth, nr;			/* Shorthand */

    max_eps = 0.0;
    nr = grids.num_r - 2;  /* Ignore silicate and ocean layers */
    nth = grids.num_theta;
    dth = grids.dth;
    dr = grids.dr;
    dth2 = dth*dth;

    /****************************************************************/

    /* Special treatment for k = 1, unequal finite difference */
    rk = grids.ice_radius[1];
    delta = dr/rk;
//    (void)fprintf(stderr,"rk %Lf delta %Lf\n",rk,delta);
    delta2 = delta*delta;
    /* Refl. BC at N pole */
    Tprev = T[1][1];
    /* Elementally */
    thi = grids.etheta[1];
    /* T[1][1] = ( T[1][2] * (2.0 + delta)
		    + T[1][0] * (4.0 - 4.0*delta)
		    + H[1][1] * dr * dr / planet.k  )
		/ (6.0 - 3.0*delta);
    */
    T[1][1] = ( T[1][2] * (1.0 + delta)
		+ T[1][0] * 2.0 * (1.0 - delta)
		+ (T[2][1] + T[1][1]) * delta2/dth2
		+ (T[2][1] - T[1][1]) * delta2 / (dth * tan(thi))
		+ H[1][1] * dr * dr / planet.k  )
	    / (3.0 - delta + 2.0*delta2/dth2);

    eps = fabs( (T[1][1] - Tprev)/Tprev );
    if (eps > max_eps) 
	max_eps = eps;

    for (i=2;i<=(nth-2);i++) {
	thi = grids.etheta[i];
	Tprev = T[i][1];
/*	T[i][1] = ( T[i][2] * (2.0 + delta)
		    + T[i][0] * (4.0 - 4.0*delta)
		    + H[i][1] * dr * dr / planet.k  )
		/ (6.0 - 3.0*delta);
*/
	T[i][1] = ( T[i][2] * (1.0 + delta)
		    + T[i][0] * 2.0 * (1.0 - delta)
		    + T[i+1][1] * delta2/dth2 * (1.0 + dth/tan(thi))
		    + T[i-1][1] * delta2/dth2 * (1.0 - dth/tan(thi))
		    + H[i][1] * dr * dr / planet.k  )
		/ (3.0 - delta + 2.0*delta2/dth2);

	eps = fabs( (T[i][1] - Tprev)/Tprev );
	if (eps > max_eps)
	    max_eps = eps;
    }    
	
    /* Refl. BC at S pole */
    Tprev = T[nth-1][1];
    /* Elementally */
    thi = grids.etheta[nth-1];
/*    T[nth-1][1] = ( T[nth-1][2] * (2.0 + delta)
		    + T[nth-1][0] * (4.0 - 4.0*delta)
		    + H[nth-1][1] * dr * dr / planet.k  )
		/ (6.0 - 3.0*delta);
*/
    T[nth-1][1] = ( T[nth-1][2] * (1.0 + delta)
		    + T[nth-1][0] * 2.0 * (1.0 - delta)
		    + (T[nth-1][1] + T[nth-2][1]) * delta2/dth2
		    + (T[nth-1][1] - T[nth-2][1]) * delta2 / (dth * tan(thi))
		    + H[nth-1][1] * dr * dr / planet.k  )
		/ (3.0 - delta + 2.0*delta2/dth2);
	    
    eps = fabs( (T[nth-1][1] - Tprev)/Tprev );
    if (eps > max_eps) 
	max_eps = eps;

    /* For consistency, make polar nodal T same as T in adjacent elem. */
    T[0][1] = T[1][1];
    T[nth][1] = T[nth-1][1];
    
    /********************************************************************/

    /* Work by layer, excluding top and bottom elements */
    for (k=2;k<=(nr-2);k++) {
	rk = grids.ice_radius[k];
	delta = dr/rk;
//	(void)fprintf(stderr,"rk %Lf delta %Lf\n",rk,delta);
	delta2 = delta*delta;
	/*
	if(k==2)
	    (void)fprintf(stderr,"%f %f %f\n %.1f %.1f %.1f %.1f %.1f \n %.1f %.1f %.1f %.1f %.1f\n", 
		0.0, delta2/dth2, 0.0,
		T[1][k],T[1][k+1],T[1][k-1],T[2][k],T[2][k],
		T[1][k+1] * (1.0 + delta2),
		T[1][k-1] * (1.0 - delta2),
		T[2][k] * delta2/dth2,
		T[2][k] * delta2/dth2,
		H[1][k] * dr * dr / properties.k);
	*/
	/* Then from N to S */
	/* Refl. BC at N pole */
	Tprev = T[1][k];
	/* Nodally */
	/* T[1][k] = ( T[1][k+1] * (1.0 + delta2)
		    + T[1][k-1] * (1.0 - delta2)
		    + 2.0 * T[2][k] * delta2/dth2
		    + H[1][k] * dr * dr / properties.k  )
		/ (2.0 + 2.0*delta2/dth2);
	*/
	/* Elementally */
	thi = grids.etheta[1];
/*	T[1][k] = (	T[1][k+1] * (1.0 + delta)
			+ T[1][k-1] * (1.0 - delta)
			+ H[1][k] * dr * dr / planet.k  )
		    / (2.0);
*/
	T[1][k] = ( T[1][k+1] * (1.0 + delta)
		    + T[1][k-1] * (1.0 - delta)
		    + (T[2][k] + T[1][k]) * delta2/dth2
		    + (T[2][k] - T[1][k]) * delta2 / (dth * tan(thi))
		    + H[1][k] * dr * dr / planet.k  )
		/ (2.0 + 2.0*delta2/dth2);

	eps = fabs( (T[1][k] - Tprev)/Tprev );
	if (eps > max_eps) 
	    max_eps = eps;

	for (i=2;i<=(nth-2);i++) {
	    thi = grids.etheta[i];
	    /*
	    if(k==2 && i<=3)
	    (void)fprintf(stderr,"%f %f %f \n %.1f %.1f %.1f %.1f %.1f \n %.1f %.1f %.1f %.1f %.1f\n", 
	        thi, delta2/dth2, dth/tan(thi),
		T[i][k],T[i][k+1],T[i][k-1],T[i+1][k],T[i-1][k],
		T[i][k+1] * (1.0 + delta2),
		T[i][k-1] * (1.0 - delta2),
		T[i+1][k] * delta2/dth2 * (1.0 + dth/tan(thi)),
		T[i-1][k] * delta2/dth2 * (1.0 - dth/tan(thi)),
		H[i][k] * dr * dr / properties.k);
	    */
	    Tprev = T[i][k];
/*	    T[i][k] = (	T[i][k+1] * (1.0 + delta)
			+ T[i][k-1] * (1.0 - delta)
			+ H[i][k] * dr * dr / planet.k  )
		    / (2.0);
*/
	    T[i][k] = (	T[i][k+1] * (1.0 + delta)
			+ T[i][k-1] * (1.0 - delta)
			+ T[i+1][k] * delta2/dth2 * (1.0 + dth/tan(thi))
			+ T[i-1][k] * delta2/dth2 * (1.0 - dth/tan(thi))
			+ H[i][k] * dr * dr / planet.k  )
		    / (2.0 + 2.0*delta2/dth2);

	    eps = fabs( (T[i][k] - Tprev)/Tprev );
	    if (eps > max_eps)
		max_eps = eps;
	}    
	
	/* Refl. BC at S pole */
	Tprev = T[nth-1][k];
	/* Nodally */
	/*
	T[nth][k] = ( T[nth][k+1] * (1.0 + delta2)
			+ T[nth][k-1] * (1.0 - delta2)
			+ 2.0 * T[nth-1][k] * delta2/dth2
			+ H[nth][k] * dr * dr / properties.k  )
		    / (2.0 + 2.0*delta2/dth2);
	 */
	/* Elementally */
	thi = grids.etheta[nth-1];
/*	T[nth-1][k] = (	T[nth-1][k+1] * (1.0 + delta)
			+ T[nth-1][k-1] * (1.0 - delta)
			+ H[nth-1][k] * dr * dr / planet.k  )
		    / (2.0);
*/
	T[nth-1][k] = ( T[nth-1][k+1] * (1.0 + delta)
		    + T[nth-1][k-1] * (1.0 - delta)
		    + (T[nth-1][k] + T[nth-2][k]) * delta2/dth2
		    + (T[nth-1][k] - T[nth-2][k]) * delta2 / (dth * tan(thi))
		    + H[nth-1][k] * dr * dr / planet.k  )
		/ (2.0 + 2.0*delta2/dth2);
	    
	eps = fabs( (T[nth-1][k] - Tprev)/Tprev );
	if (eps > max_eps) 
	    max_eps = eps;

	/* For consistency, make polar nodal T same as T in adjacent elem. */
	T[0][k] = T[1][k];
	T[nth][k] = T[nth-1][k];
    }

    /********************************************************************/

    /* Special treatment for k = nr-1, unequal finite difference */
    rk = grids.ice_radius[nr-1];
    delta = dr/rk;
//    (void)fprintf(stderr,"rk %Lf delta %Lf\n",rk,delta);
    delta2 = delta*delta;
    /* Refl. BC at N pole */
    Tprev = T[1][nr-1];
    /* Elementally */
    thi = grids.etheta[1];
/*    T[1][nr-1] = ( T[1][nr] * (4.0 + 4.0*delta)
		+ T[1][nr-2] * (2.0 - delta)
		+ H[1][nr-1] * dr * dr / planet.k  )
	    / (6.0 + 3.0*delta);
*/
    T[1][nr-1] = ( T[1][nr] * 2.0 * (1.0 + delta)
		+ T[1][nr-2] * (1.0 - delta)
		+ (T[2][nr-1] + T[1][nr-1]) * delta2/dth2
		+ (T[2][nr-1] - T[1][nr-1]) * delta2 / (dth * tan(thi))
		+ H[1][nr-1] * dr * dr / planet.k  )
	    / (3.0 + delta + 2.0*delta2/dth2);

    eps = fabs( (T[1][nr-1] - Tprev)/Tprev );
    if (eps > max_eps) 
	max_eps = eps;

    for (i=2;i<=(nth-2);i++) {
	thi = grids.theta[i];
	Tprev = T[i][nr-1];
/*    	T[i][nr-1] = ( T[i][nr] * (4.0 + 4.0*delta)
		+ T[i][nr-2] * (2.0 - delta)
		+ H[i][nr-1] * dr * dr / planet.k  )
	    / (6.0 + 3.0*delta);
*/
	T[i][nr-1] = ( T[i][nr] * 2.0 * (1.0 + delta)
		    + T[i][nr-2] * (1.0 - delta)
		    + T[i+1][nr-1] * delta2/dth2 * (1.0 + dth/tan(thi))
		    + T[i-1][nr-1] * delta2/dth2 * (1.0 - dth/tan(thi))
		    + H[i][nr-1] * dr * dr / planet.k  )
		/ (3.0 + delta + 2.0*delta2/dth2);

	eps = fabs( (T[i][nr-1] - Tprev)/Tprev );
	if (eps > max_eps)
	    max_eps = eps;
    }    
	
    /* Refl. BC at S pole */
    Tprev = T[nth-1][nr-1];
    /* Elementally */
    thi = grids.etheta[nth-1];
/*    T[nth-1][nr-1] = ( T[nth-1][nr] * (4.0 + 4.0*delta)
		+ T[nth-1][nr-2] * (2.0 - delta)
		+ H[nth-1][nr-1] * dr * dr / planet.k  )
	    / (6.0 + 3.0*delta);
*/
    T[nth-1][nr-1] = ( T[nth-1][nr] * 2.0 * (1.0 + delta)
		    + T[nth-1][nr-2] * (1.0 - delta)
		    + (T[nth-1][nr-1] + T[nth-2][nr-1]) * delta2/dth2
		    + (T[nth-1][nr-1] - T[nth-2][nr-1]) * delta2 / (dth * tan(thi))
		    + H[nth-1][nr-1] * dr * dr / planet.k  )
		/ (3.0 + delta + 2.0*delta2/dth2);
	    
    eps = fabs( (T[nth-1][nr-1] - Tprev)/Tprev );
    if (eps > max_eps) 
	max_eps = eps;

    /* For consistency, make polar nodal T same as T in adjacent elem. */
    T[0][nr-1] = T[1][nr-1];
    T[nth][nr-1] = T[nth-1][nr-1];
    
    return(max_eps);

}
