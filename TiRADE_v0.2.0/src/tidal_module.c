/********************************************************
 * tidal_module -- Routines to compute the spatially-	*
 *	varying tidal heating within a spherically			*
 *	symmetric planet or moon composed of an 			*
 *	arbitrary number of layers.							*
 *														*
 * Purpose:												*
 * 	Originally intended as input parameter for a		*
 *	convection program such as Citcom (i.e. 			*
 *	spatially varying internal heating, but the			*
 *	average properties can be analyzed on their 		*
 *	own.												*
 *														*
 * Notes:												*
 *	Modified from a free-standing code,					*
 *	tidal_heating										*
 *														*
 * Author:  James Roberts								*
 *														*
 * Last modified:  31 Jan 2014							*
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#ifndef __TIRADE_H__
#define __TIRADE_H__
#include "tirade.h"
#endif

#include "utility.h"

#define L 2			/* spherical harmonic degree */
#define NUM_T 8		/* Total number of timesteps for heating average */

#define SWAP(a_elem,b_elem) {temp=(a_elem); (a_elem)=(b_elem); (b_elem)=temp;}

extern long double **heating_axi_ice;
extern struct grid_info grids;
extern struct globals planet;
extern struct matprop *layers;

/************************************************************
 * tidal_heating -- Main procedure and Wrapper for all the  *
 *	routines to compute the tidal heating.		            *
 *							                                *
 * Parameters						                        *
 *	max_layers -- number of layers in the planet	        *
 *	output_flags -- controls which data are output	        *
 *	directory -- output directory			                *
 ************************************************************/
void tidal_heating(int max_layers, char output_flags, char directory[], 
	int tc_iter, char* true_fluid_flag)
{
    int i,j,k,ti;					/* counters */
    complex long double *y_i[7];	/* Radial functions from S&V 2004 */
    long double **potential;		/* Tidal Potential from Kaula 1964 */
    long double **dpot[5];			/* and its derivatives */
    long double time;
    complex long double Loveh2, Lovek2;	/* Love Number */
    long double phaseh2, phasek2;		/* Phase Lags */
    complex long double k2an;
    long double Ediss;					/* Total Dissipated Power */
    long double Q, Qz;					/* Quality Factor */
    long double maxwell_time;
    complex long double *strain[7];		/* Strain Tensor */
    complex long double *stress[7];		/* Stress Tensor */
    long double *tidal_heat;			/* Tidal heating rate */
    long double SA;
    char filename[80];

    FILE *out_file1;
    FILE *out_file2;

    void radial_functions(int, complex long double *y_i[], char*);
    void tidal_potential(long double, long double **, long double **dpot[]);
    void stress_strain(int, complex long double *y_i[], long double **, 
	    long double **dpot[], complex long double *strain[], 
	    complex long double *stress[]);
    void heating(complex long double *strain[], complex long double *stress[],
	    long double *);
    void output_data(int, complex long double *y_i[], long double **, 
	    long double **dpot[], complex long double *strain[],
	    complex long double *stress[], long double *tidal_heat, char,
	    char directory[], int);

    /* Allocate Memory */

    for(i=1;i<=6;i++) 
    	y_i[i] = (complex long double *)malloc((grids.num_r+1)
    			*sizeof(complex long double));

    potential = (long double **)malloc((grids.num_theta+1)*sizeof(long double));
    potential[0] = (long double *)malloc((grids.num_theta+1)*(grids.num_phi+1)
	    *sizeof(long double));
    	for(i=1;i<=grids.num_theta;i++)
    		potential[i] = potential[i-1] + (grids.num_phi+1);

    for(j=0;j<=4;j++) {
    	dpot[j] = (long double **)malloc((grids.num_theta+1)*sizeof(long double));
    	dpot[j][0] = (long double *)malloc((grids.num_theta+1)*(grids.num_phi+1)
    			*sizeof(long double));
    	for(i=1;i<=grids.num_theta;i++)
    		dpot[j][i] = dpot[j][i-1] + (grids.num_phi+1);
    }

    for(i=1; i<=6; i++) {
    	strain[i] = (complex long double *)malloc((grids.num_theta+1)
    			*(grids.num_phi+1)*(grids.num_r+1)*sizeof(complex long double));
    	stress[i] = (complex long double *)malloc((grids.num_theta+1)
    			*(grids.num_phi+1)*(grids.num_r+1)*sizeof(complex long double));
    }

    tidal_heat = (long double *)malloc((grids.num_theta+1)*(grids.num_phi+1)
	    *(grids.num_r+1)*sizeof(long double));

    /* Verify new viscosities */
    /*for (k=0;k<=(max_layers-1);k++)
	(void)fprintf(stderr,"layer %d visc %Le Pa s\n",k,layers[k].visc);
*/
    /* Get the radial functions y_i at each depth ... */

    radial_functions(max_layers,y_i,true_fluid_flag);

    /* Find the Love Number and total dissipated energy*/

    k2an = 1.5 / (1.0 + 19.0*layers[max_layers-1].rigidity / (2.0 * layers[max_layers-1].dens * layers[max_layers-1].g * layers[max_layers-1].rad)); 
    (void)fprintf(stdout,"k2an = %Lg + %Lgi\n",creall(k2an),cimagl(k2an));

    Lovek2 = (-y_i[5][grids.num_r] - 1.0);
    Loveh2 = y_i[1][grids.num_r] * layers[max_layers-1].g;
    phaseh2 = atan(-cimagl(Loveh2)/creall(Loveh2)) * 180.0 / PI;
    phasek2 = atan(-cimagl(Lovek2)/creall(Lovek2)) *180.0 / PI;
    //(void)fprintf(stdout,"g = %Lg m s^-2\n",layers[max_layers-1].g);
    (void)fprintf(stdout,"h2 = %Lg + %Lg, amplitude %Lg, phase %Lg\n",
	    creall(Loveh2),cimagl(Loveh2),cabsl(Loveh2),phaseh2);
    (void)fprintf(stdout,"k2 = %Lg + %Lg, amplitude %Lg, phase %Lg\n",
	    creall(Lovek2),cimagl(Lovek2),cabsl(Lovek2),phasek2);
/*    (void)fprintf(stdout,"k2/h2 = %Lg + %Lgi, phase = %Lg %Lg\n",creall(Lovek2/Loveh2),cimagl(Lovek2/Loveh2),phaseh2,phasek2);
*/
    Ediss = -10.5*cimagl(Lovek2) 
	    * pow((planet.freq * layers[max_layers-1].rad),5)
	    * planet.ecc * planet.ecc / G;
    SA = 4.0*PI* layers[max_layers-1].rad*layers[max_layers-1].rad;
    maxwell_time = layers[max_layers-1].visc/layers[max_layers-1].shear;
    Q = (1.0 + planet.freq*maxwell_time * planet.freq*maxwell_time) 
	/ (2.0 * planet.freq*maxwell_time);

    (void)fprintf(stdout,"Ediss = %Lg W, ",Ediss);
    (void)fprintf(stdout,"Ediss/Vol = %Lg W/m^3\n",Ediss/planet.vol);
    //(void)fprintf(stdout,"HF = %Lg W/m^2\n",Ediss/SA);
    (void)fprintf(stdout,"mu = %Lg + %Lgi\n",creall(layers[max_layers-1].rigidity),cimagl(layers[max_layers-1].rigidity));
    //(void)fprintf(stdout,"tau = %Lg s\n",maxwell_time);
    //(void)fprintf(stdout,"Q = %Lg\n",Q);
/*    out_file1 = fopen("MS/MSfig3.dat","a");
    (void)fprintf(out_file1,"%Lg %Lg %Lg\n",(layers[max_layers-1].rad - layers[max_layers-3].rad)/1000.0, layers[max_layers-2].visc, phaseh2);
    (void)fclose(out_file1);
*/
    /* 
     * We care about the heating averaged over one forcing cycle.
     * Therefore, we need to do the following steps at a number of
     * timesteps in order to get an average.
     */

    (void)fprintf(stderr,"ti =  ");
    for(ti=0; ti<NUM_T; ti++) {
    	time = ti/(long double)NUM_T * (2.0*PI)/planet.freq;
    	(void)fprintf(stderr,"%d, ",ti);

    	/* Get Tidal Potential; lateral variation */
    
    	tidal_potential(time,potential,dpot);

    	/* Find stress and strain tensors, combine lateral and radial
	   	   variations */

    	stress_strain(max_layers,y_i,potential,dpot,strain,stress);

    	/* Find heating */

    	heating(strain,stress,tidal_heat);
    }
	(void)fprintf(stderr,"\n");

    /* Write the heating to a file for use by others */
    /* Pass in other data, why not? */
    output_data(max_layers,y_i,potential,dpot,strain,stress,tidal_heat,
	    output_flags, directory, tc_iter);

    /* Free Memory */

    for(i=1; i<=6; i++) 
    	(void)free((complex long double *)y_i[i]);

    for(i=1; i<=6; i++) {
    	(void)free((complex long double *)strain[i]);
    	(void)free((complex long double *)stress[i]);
    }

    (void)free((long double *)potential[0]);
    (void)free((long double **)potential);

    for(j=0;j<=4;j++) {
    	(void)free((long double *)dpot[j][0]);
    	(void)free((long double **)dpot[j]);
    }
    
    (void)free((long double *)tidal_heat);

}


/************************************************************
 * radial_functions -- given the material properties for    *
 *		each layer, applies a propogator matrix	            *
 *		method (Sabadini and Vermeersen, 2004) to           *
 *		solve for the radial dependence of	                *
 *		displacement, stress, potential, and	            *
 *		potential gradient.			                        *
 *		We assume incompressibility.		                *
 *							                                *
 * Parameters			                                    *
 *	max_layers -- number of layers in the planet	        *
 *	y_i -- matrix containing the solution for the	        *
 *		radial functions. This is empty when passed         *
 *		in.					                                *
 ************************************************************/

void radial_functions(int max_layers, complex long double *y_i[], char* true_fluid_flag)
{
    long double *Ic[7];		    		/* Interface Matrix at CMB */
    complex long double *Y[7][7];	    /* Propagator Matrix */
    complex long double *Yinv[7][7]; 
    complex long double *Ylayer[7];
    complex long double *Yinvlayer[7];
    complex long double *propmat[7]; 
    complex long double *product[7]; 
    int i,j,k;			    			/* matrix counters */
    int layer;			    			/* layer counter */
    int current_layer;
    long double r,rho,g;	    	/* = layers[layer].vars, for convenience */
    complex long double mu;
    complex long double Cc[4];	    	/* Conditions at top of "fluid" core */

    void get_layer_propmat(int, long double, complex long double *Ylayer[], 
	    complex long double *Yinvlayer[]);
	void get_true_fluid_propmat(int, long double, complex long double *Ylayer[], complex long double *Yinvlayer[]);
    void mat_multi_mat(complex long double *Ylayer[], 
	    complex long double *propmat[], complex long double *product[], 
	    int, int, int);
    void apply_boundary_conditions(complex long double *Y[][7], 
	    complex long double *Yinv[][7], long double *Ic[4], 
	    complex long double Cc[],int);
    int find_layer(long double, int);
    
    /* Allocate Memory */
    for (i=1; i<=6; i++) {
    	Ic[i] = (long double *)malloc(4*sizeof(long double));
    	Yinvlayer[i] = (complex long double *)malloc(7*sizeof(complex long double));
    	Ylayer[i] = (complex long double *)malloc(7*sizeof(complex long double));
    	propmat[i] = (complex long double *)malloc(7*sizeof(complex long double));
    	product[i] = (complex long double *)malloc(7*sizeof(complex long double));

    	for (j=1; j<=6; j++) {
    		Y[i][j] = (complex long double *)malloc((max_layers+2)*sizeof(complex long double));
    		Yinv[i][j] = (complex long double *)malloc((max_layers+2)*sizeof(complex long double));
    	}

    }

    /* Set up Ic, interface matrix at CMB */

    Ic[1][1] = -pow(layers[0].rad,(L-1)) / layers[0].A;
    Ic[1][2] = 0.0;
    Ic[1][3] = 1.0;

    Ic[2][1] = 0.0;
    Ic[2][2] = 1.0;
    Ic[2][3] = 0.0;

    Ic[3][1] = 0.0;
    Ic[3][2] = 0.0;
    Ic[3][3] = layers[0].dens * layers[0].A * layers[0].rad;

    Ic[4][1] = 0.0;
    Ic[4][2] = 0.0;
    Ic[4][3] = 0.0;

    Ic[5][1] = pow(layers[0].rad,L);
    Ic[5][2] = 0.0;
    Ic[5][3] = 0.0;

    Ic[6][1] = 2*(L-1) * pow(layers[0].rad,(L-1));
    Ic[6][2] = 0.0;
    Ic[6][3] = 3.0 * layers[0].A;

    /*
     * Assemble Y matrix and its inverse for each layer.
     * Careful, L is an integer. Need denom. to be long double for division
     * to work right.
     */
    for (layer=1; layer<max_layers; layer++) {
		r = layers[layer].rad;
		if (*true_fluid_flag && creall(layers[layer].visc) < planet.true_fluid_cutoff_visc && creall(layers[layer].shear) < planet.true_fluid_cutoff_rigidity && layer > 1) { 
			get_true_fluid_propmat(layer, r, Ylayer, Yinvlayer);
		}
		else { get_layer_propmat(layer,r,Ylayer,Yinvlayer); }
		for (i = 1; i <= 6; i++) {
			for (j = 1; j <= 6; j++) {
				Y[i][j][layer] = Ylayer[i][j];
				Yinv[i][j][layer] = Yinvlayer[i][j];
			}
		}
    }

    /* 
     * Apply the Boundary Conditions 
     * 
     * Ic describes the conditions at the CMB, assuming a fluid core
     * (which can be arbitrily small), in terms of a set of constants
     * Cc, related to the conditions at the surface.  
     * We use the propagator matrix P = Product(Y*Y') to get
     * the conditions in the top layer in terms of Cc,
     * then directly apply the surface boundary conditions to solve 
     * for Cc.  We can then use Cc to solve for y_i at any point .
     */

    apply_boundary_conditions(Y,Yinv,Ic,Cc,max_layers);
    /*for(i=1;i<=3;i++)
	(void)fprintf(stderr,"Cc[%d] = %Le + %Le i\n",
		i,creall(Cc[i]),cimagl(Cc[i]));
    */
    /*
     * Solve for radial functions, y_i  
     *
     * Now we have the vector Cc. We can solve for y_i at any
     * point r. We use the pre-existing Y, Yinv matrices up to
     * the layer below the layer in which point r resides. We then
     * split that layer at radius r and Propagate up to that point.
     */

    for(k=1;k<=grids.num_r;k++) {
    	/* Calc yi at center of each layer, and also at the surface */
    	r = grids.eradius[k];	    /* Interior Element */

    	/* Find which layer point r is in */
    	current_layer = find_layer(r,max_layers);

    	for(i=1;i<=6;i++) {
    		for(j=1;j<=3;j++) {
    			propmat[i][j] = Ic[i][j];
    		}
    	}

    	/* Get Prop Mat up to current layer */
    	for(layer=1; layer<=(current_layer-1); layer++) {

    		/* Get Y, Y' for layer */
    		for(i=1;i<=6;i++) {
    			for(j=1;j<=6;j++) {
    				Ylayer[i][j] = Y[i][j][layer];
    				Yinvlayer[i][j] = Yinv[i][j][layer];
    			}
    		}
	
    		/* product = Y'(layer) * propmat */
    		mat_multi_mat(Yinvlayer,propmat,product,6,3,6);
    		for(i=1;i<=6;i++)
    			for(j=1;j<=3;j++)
    				propmat[i][j] = product[i][j];
	
    		/* product = Y(layer) * propmat */
    		mat_multi_mat(Ylayer,propmat,product,6,3,6);
    		for(i=1;i<=6;i++)
    			for(j=1;j<=3;j++)
    				propmat[i][j] = product[i][j];
    	}

    	/* Add in Partial Layer */
    	//get_layer_propmat(layer,r,Ylayer,Yinvlayer);
		if (*true_fluid_flag && creall(layers[layer].visc) < planet.true_fluid_cutoff_visc && creall(layers[layer].shear) < planet.true_fluid_cutoff_rigidity && layer > 1) {
			get_true_fluid_propmat(layer, r, Ylayer, Yinvlayer);
		}
		else { get_layer_propmat(layer,r,Ylayer,Yinvlayer); }

    	/* product = Y'(layer) * propmat */
    	mat_multi_mat(Yinvlayer,propmat,product,6,3,6);

    	for(i=1;i<=6;i++) {
    		for(j=1;j<=3;j++) {
    			propmat[i][j] = product[i][j];
    		}
    	}
	
    	/* product = Y(layer) * propmat */
    	mat_multi_mat(Ylayer,propmat,product,6,3,6);

    	for(i=1;i<=6;i++) {
    		for(j=1;j<=3;j++) {
    			propmat[i][j] = product[i][j];
    		}
    	}

    	/* find y_i(r) = P * Cc */
    	/*
    	 * Note that the Im components of y3 and particularly y4 may	*
    	 * disagree with the MM solution by as much as 20%.  This is	*
    	 * probably due to the low magnitude of these components	*
    	 * relative to the reals.  Only the magnitude is off; the shape	*
    	 * agrees quite well.  The real components of y3 and y4, and	*
    	 * all components of the other yi's have excellent agreement.	*
    	 */
    	for(i=1;i<=6;i++) {
    		y_i[i][k] = 0.0;	/* Initialize y_i */
    		for(j=1;j<=3;j++)
    			y_i[i][k] += propmat[i][j] * Cc[j];
    	}

    }

    /* Free Memory */
    for (i=1; i<=6; i++) {
    	free((long double *)Ic[i]);
    	free((complex long double *)Yinvlayer[i]);
    	free((complex long double *)Ylayer[i]);
    	free((complex long double *)propmat[i]);
    	free((complex long double *)product[i]);

    	for (j=1; j<=6; j++) {
    		free((complex long double *)Y[i][j]);
    		free((complex long double *)Yinv[i][j]);
    	}
    }

}

void get_true_fluid_propmat(int layer, long double r, complex long double *Ylayer[], complex long double *Yinvlayer[]) {
	// mu = 0
	// propagation of mechanical properties, y1, y3, y4 and potential stress y6 basically go away
	// only left w/ y5, y7
	void diag_multi_mat(complex long double *Yinvlayer[], long double D[], complex long double *Ybar[], int);

	// note to self SEE PAGE 61 (Sabadini)
	int i, j;
	/*for (i = 1; i < 7; i++) {
		for (j = 1; j < 7; j++) {
			Ylayer[i][j] = 0.0;
		}
	}*/

	// d(y5)/dr = y7
	// d(y7)/dr = (l(l+1)/r^2)y5 - (2/r)y7
	long double rho = layers[layer].dens;
    long double g = layers[layer].g;

	// NOTE: B = Y(r_n, s)(Y(r_{n + 1}, s))^-1
	// NOTE: Ylayer[row][col]

	// using radius at ocean-ice interface

	// y1
	Ylayer[1][1] = L * pow(r, L + 1) / (long double) (4 * L + 6);
	Ylayer[1][2] = pow(r, L - 1);
	Ylayer[1][3] = Ylayer[1][6] = 0.0;
	Ylayer[1][4] = (L + 1) * pow(r, -L) / (long double) (4 * L - 2);
	Ylayer[1][5] = pow(r, -L - 2);

	// y2
	Ylayer[2][1] = (L + 3) * pow(r, L + 1) / (long double) (2 * (2 * L + 3) * (L + 1));
	Ylayer[2][2] = pow(r, L - 1) / (long double) L;
	Ylayer[2][4] = (2 - L) * pow(r, -L) / (long double) (2 * L * (2 * L - 1));
	Ylayer[2][5] = -pow(r, -L - 2)/ (long double) (L + 1);

	// y3
	Ylayer[3][1] = (L * rho * g * r) * pow(r, L) / (long double) (4 * L + 6);
	Ylayer[3][2] = (rho * g * r) * pow(r, L-2);
	Ylayer[3][3] = rho * pow(r, L);
	Ylayer[3][4] = (L + 1) * rho * g * r / (long double) (2 * (2 * L - 1) * pow(r, L + 1));
	Ylayer[3][5] = (rho * g * r) / (long double) pow(r, L + 3);
	Ylayer[3][6] = rho / (long double) pow(r, L + 1);

	// y5
	Ylayer[5][3] = pow(r, L);
	Ylayer[5][6] = 1 / (long double) pow(r, L + 1);

	// y6
	Ylayer[6][1] = 2 * PI * G * rho * L * pow(r, L + 1) / (long double) (2 * L + 3);
	Ylayer[6][2] = 4 * PI * G * rho * pow(r, L - 1);
	Ylayer[6][3] = (2 * L + 1) * pow(r, L - 1);
	Ylayer[6][4] = (2 * PI * G * rho * (L + 1)) / (long double) ((2 * L - 1) * pow(r, L));
	Ylayer[6][5] = (4 * PI * G * rho) / (long double) pow(r, L + 2);

	// Now, define Yinv... we need D and Ybar for this
	r = layers[layer - 1].rad; // get radius at mantle-ocean interface

	// note to self - don't change D[i]
	long double D[7];
	D[1] = (long double)(L+1) / (long double)pow(r,(L+1));
    D[2] = (long double)(L*(L+1)) / (long double)(2*(2*L-1)*pow(r,(L-1)));
    D[3] = 1.0 / (long double)pow(r,(L-1));
    D[4] = (long double)(L*pow(r,L));
    D[5] = ( (long double)(L*(L+1)) / (long double)(2*(2*L+3)) ) * pow(r,(L+2));
    D[6] = -(long double)(pow(r,(L+1)));

	for(i=1;i<=6;i++) D[i] /= (long double)(2*L+1);

	// Allocate memory for Ybar
	complex long double *Ybar[7];
	for (i = 1; i <= 6; i++) Ybar[i] = (complex long double *) malloc(7 * sizeof(complex long double));

	// Assign Ybar
	long double complex mu = layers[layer - 1].rigidity;
	//rho = layers[layer - 1].dens;
	//g = layers[layer - 1].g;

	Ybar[1][1] = rho*g*r/mu - 2*(L+2);
	Ybar[1][2] = (long double)(2*L*(L+2));
    Ybar[1][3] = -r/mu;
    Ybar[1][4] = L*r/mu;
    Ybar[1][5] = rho*r/mu;
    Ybar[1][6] = 0.0;
    
    Ybar[2][1] = -rho*g*r/mu + (long double)(2*(L*L + 3*L - 1)) / (long double)(L+1);
    Ybar[2][2] = (long double)(-2*(L*L-1));
    Ybar[2][3] = r/mu;
    Ybar[2][4] = ((2-L)*r)/mu;
    Ybar[2][5] = -rho*r/mu;
    Ybar[2][6] = 0.0;

    Ybar[3][1] = 4.0*PI*G*rho;
    Ybar[3][2] = 0.0;
    Ybar[3][3] = 0.0;
    Ybar[3][4] = 0.0;
    Ybar[3][5] = 0.0;
    Ybar[3][6] = -1.0;

    Ybar[4][1] = rho*g*r/mu + (long double)(2*(L-1)) ;
    Ybar[4][2] = (long double)(2*(L*L-1));
    Ybar[4][3] = -r/mu;
    Ybar[4][4] = -((L+1)*r)/mu;
    Ybar[4][5] = rho*r/mu;
    Ybar[4][6] = 0.0;
    
    Ybar[5][1] = -rho*g*r/mu - (long double)(2*(L*L - L - 3)) / (long double)(L);
    Ybar[5][2] = (long double)(-2*L*(L+2));
    Ybar[5][3] = r/mu;
    Ybar[5][4] = ((L+3)*r)/mu;
    Ybar[5][5] = -rho*r/mu;
    Ybar[5][6] = 0.0;
	
    Ybar[6][1] = 4.0*PI*G*rho*r;
    Ybar[6][2] = 0.0;
    Ybar[6][3] = 0.0;
    Ybar[6][4] = 0.0;
    Ybar[6][5] = (long double)(2*L+1);
    Ybar[6][6] = -r;

	// Define Yinvlayer
	diag_multi_mat(Yinvlayer,D,Ybar,6);

	// Deallocate heap memory
	for (i = 1; i <= 6; i++) free(Ybar[i]);
	return;
}

/************************************************************
 * get_layer_propmat -- a routine to acquire the Y and Yinv *
 *	    matrices for a layer given the physical	            *
 *	    properties of the layer			                    *
 *							                                *
 * Parameters						                        *
 *	    layer   --  index of the layer		                *
 *	    Y	    --  Y for the layer at r		            *
 *	    Yinv    --	Y' for the layer at r_inner	            *
 ************************************************************/

void get_layer_propmat(int layer, long double r, complex long double *Ylayer[], 
			complex long double *Yinvlayer[])
{
    complex long double *Ybar[7];
    long double D[7];
    long double rho, g;	/* placeholder variables */
    complex long double mu;
    int i;		/* Matrix counter */

    void diag_multi_mat(complex long double *Yinvlayer[], long double D[], 
	    complex long double *Ybar[], int);

    for(i=1;i<=6;i++)
    	Ybar[i] = (complex long double *)malloc(7*sizeof(complex long double));

    rho = layers[layer].dens;
    mu = layers[layer].rigidity;
    g = layers[layer].g;

    /* Y matrix */

    Ylayer[1][1] = (L*pow(r,(L+1))) / (2.0*(2*L+3));
    Ylayer[1][2] = pow(r,(L-1));
    Ylayer[1][3] = 0.0;
    Ylayer[1][4] = ((L+1)*pow(r,(-L))) / (2.0*(2*L-1));
    Ylayer[1][5] = pow(r,(-L-2));
    Ylayer[1][6] = 0.0;

    Ylayer[2][1] = ((L+3)*pow(r,(L+1))) / (2.0*(2*L+3)*(L+1));
    Ylayer[2][2] = pow(r,(L-1))/L;
    Ylayer[2][3] = 0.0;
    Ylayer[2][4] = ((2-L)*pow(r,(-L))) / (2.0*L*(2*L-1));
    Ylayer[2][5] = -1.0 * pow(r,(-L-2)) / (long double)(L+1);
    Ylayer[2][6] = 0.0;

    Ylayer[3][1] = ((L*rho*g*r + 2*(L*L-L-3)*mu) * pow(r,L)) /
		    (long double)(2*(2*L+3));
    Ylayer[3][2] = (rho*g*r + 2*(L-1)*mu) * pow(r,(L-2));
    Ylayer[3][3] = -rho * pow(r,L);
    Ylayer[3][4] = ((L+1)*rho*g*r - 2*(L*L+3*L-1)*mu) /
    		    (long double)(2*(2*L-1)*pow(r,(L+1)));
    Ylayer[3][5] = (rho*g*r - 2*(L+2)*mu) / (long double)pow(r,(L+3));
    Ylayer[3][6] = -1.0 * rho / pow(r,(L+1));
	
    Ylayer[4][1] = (L*(L+2)*mu*pow(r,L)) / (long double)((2*L+3)*(L+1));
    Ylayer[4][2] = (2*(L-1)*mu*pow(r,(L-2))) / (long double)L;
    Ylayer[4][3] = 0.0;
    Ylayer[4][4] = ((L*L-1)*mu) / (long double)(L*(2*L-1)*pow(r,(L+1)));
    Ylayer[4][5] = (2*(L+2)*mu) / (long double)((L+1)*pow(r,(L+3)));
    Ylayer[4][6] = 0.0;
    
    Ylayer[5][1] = 0.0;
    Ylayer[5][2] = 0.0;
    Ylayer[5][3] = -1.0 * pow(r,L);
    Ylayer[5][4] = 0.0;
    Ylayer[5][5] = 0.0;
    Ylayer[5][6] = -1.0 / (long double)pow(r,(L+1));
	
    Ylayer[6][1] = (2.0*PI*G*rho*L*pow(r,(L+1))) / (long double)(2*L+3);
    Ylayer[6][2] = 4.0*PI*G*rho*pow(r,(L-1));
    Ylayer[6][3] = -1.0 * (2*L+1)*pow(r,(L-1));
    Ylayer[6][4] = (2.0*PI*G*rho*(L+1)) / (long double)((2*L-1)*pow(r,L));
    Ylayer[6][5] = 4.0*PI*G*rho / (long double)pow(r,(L+2));
    Ylayer[6][6] = 0.0;

    /* 
     * Get Y'.  This is wanted at the bottom of each layer, not the top,
     * so use the outer radius of the layer below.
     * For a benchmark test YY' = I, the same radius must be used.
     */

    r = layers[layer-1].rad;  /* Inner radius of the layer */

    /* Diagonal Matrix D (used for inversion of Y) */

    D[1] = (long double)(L+1) / (long double)pow(r,(L+1));
    D[2] = (long double)(L*(L+1)) / (long double)(2*(2*L-1)*pow(r,(L-1)));
    D[3] = 1.0 / (long double)pow(r,(L-1));
    D[4] = (long double)(L*pow(r,L));
    D[5] = ( (long double)(L*(L+1)) / (long double)(2*(2*L+3)) ) * pow(r,(L+2));
    D[6] = -(long double)(pow(r,(L+1)));

    for(i=1;i<=6;i++)
        D[i] /= (long double)(2*L+1);

    /* Ybar (used for inversion of Y) */
    Ybar[1][1] = rho*g*r/mu - (long double)(2*(L+2));
    Ybar[1][2] = (long double)(2*L*(L+2));
    Ybar[1][3] = -r/mu;
    Ybar[1][4] = L*r/mu;
    Ybar[1][5] = rho*r/mu;
    Ybar[1][6] = 0.0;
    
    Ybar[2][1] = -rho*g*r/mu + (long double)(2*(L*L + 3*L - 1)) / (long double)(L+1);
    Ybar[2][2] = (long double)(-2*(L*L-1));
    Ybar[2][3] = r/mu;
    Ybar[2][4] = ((2-L)*r)/mu;
    Ybar[2][5] = -rho*r/mu;
    Ybar[2][6] = 0.0;

    Ybar[3][1] = 4.0*PI*G*rho;
    Ybar[3][2] = 0.0;
    Ybar[3][3] = 0.0;
    Ybar[3][4] = 0.0;
    Ybar[3][5] = 0.0;
    Ybar[3][6] = -1.0;

    Ybar[4][1] = rho*g*r/mu + (long double)(2*(L-1)) ;
    Ybar[4][2] = (long double)(2*(L*L-1));
    Ybar[4][3] = -r/mu;
    Ybar[4][4] = -((L+1)*r)/mu;
    Ybar[4][5] = rho*r/mu;
    Ybar[4][6] = 0.0;
    
    Ybar[5][1] = -rho*g*r/mu - (long double)(2*(L*L - L - 3)) / (long double)(L);
    Ybar[5][2] = (long double)(-2*L*(L+2));
    Ybar[5][3] = r/mu;
    Ybar[5][4] = ((L+3)*r)/mu;
    Ybar[5][5] = -rho*r/mu;
    Ybar[5][6] = 0.0;
	
    Ybar[6][1] = 4.0*PI*G*rho*r;
    Ybar[6][2] = 0.0;
    Ybar[6][3] = 0.0;
    Ybar[6][4] = 0.0;
    Ybar[6][5] = (long double)(2*L+1);
    Ybar[6][6] = -r;

    /* Get Yinv = D * Ybar; */
    diag_multi_mat(Yinvlayer,D,Ybar,6);

    for(i=1;i<=6;i++)
    	free((complex long double *)Ybar[i]);
}

/************************************************************
 * diag_multi_mat-- a routine to multiply a diagonal matrix *
 *		by a square matrix.		                            *
 *					                                        *
 *		A = D * B				                            *
 *					                                        *
 * Parameters				                                *
 *	a_mat -- n by n square matrix to contain the result     *
 *	d_vec -- the diagonal of a square diagonal matrix.      *
 *	    This is passed as a vector of length n.	            *
 *	b_mat -- A square matrix			                    *
 *	n -- the size of the matrices			                *
 ************************************************************/

void diag_multi_mat(complex long double *a_mat[], long double d_vec[],
	complex long double *b_mat[], int n)
{
    int i,j;	/* Matrix indices */

    for (i=1; i<=n; i++)
	for (j=1; j<=n; j++)
	    a_mat[i][j] = d_vec[i] * b_mat[i][j];
}

/************************************************************
 * mat_multi_mat -- a routine to multiply two matrices.		*
 * 		This is presently too dumb to multiply a matrix by 	*
 * 		vector. Each matrix must be a 2D array.				*
 *							                                *
 *		C = A * B				                            *
 *							                                *
 * Parameters						                        *
 *	a_mat -- a matrix (m x o)		     		            *
 *	b_mat -- a matrix (o x n)				                *
 *	c_mat -- the product (m x n)			                *
 *	m -- rows in a_mat				                        *
 *	n -- cols in b_mat				                        *
 *	o -- cols in a_mat = rows in b_mat		                *
 ************************************************************/

void mat_multi_mat(complex long double *a_mat[], complex long double *b_mat[], complex long double *c_mat[], int m, int n, int o)
{
    int i,j;	/* Matrix indices */
    int a;	/* Input matrix index */

    /* Initialize C to 0 */
    for (i=1; i<=m; i++)
    	for (j=1; j<=n; j++)
    		c_mat[i][j] = 0.0;

    /* C = A * B */
    for (i=1; i<=m; i++)
    	for (j=1; j<=n; j++)
    		for(a=1; a<=o; a++) {
    			c_mat[i][j] += a_mat[i][a] * b_mat[a][j];
    		}

}

/************************************************************
 * apply_boundary_conditions -- Conditions at the CMB are   *
 *		known in terms of a vector Cc which depends         *
 *		on the surface BCs.  This function uses the         *
 *		Prop. Mat. method to solve for Cc, so that          *
 *		we can get the actual conditions for each           *
 *		layer.					                            *
 *							                                *
 *		P_k = Prod[Y(r_layer) Y'(r_(layer-1))] * Ic         *
 *							                                *
 *		y_i = P_k * Cc				                        *
 *							                                *
 * Parameters						                        *
 *		Y	-- Fundamental Matrix		                    *
 *		Yinv	-- Inverse of Y			                    *
 *		Cc	-- Conditions at CMB		                    *
 ************************************************************/

void apply_boundary_conditions(complex long double *Y[][7], 
	complex long double *Yinv[][7], long double *Ic[], 
	complex long double Cc[], int max_layers)
{
    complex long double *Ylayer[7];	    /* Y matrix at top of layer */
    complex long double *Yinvlayer[7];	    /* Yinv at bottom of layer */
    complex long double *product[7];	    /* Result of matrix mult */
    complex long double *propmat[7];	    /* Propagator matrix */
    complex long double *surfs[7];	    
    long double BCs[4];			    /* Surface Boundary Conditions */
    int i,j;				    /* Matrix Counters */
    int layer;				    /* Layer Counter */

    void mat_multi_mat(complex long double *Ylayer[], 
	    complex long double *Yinvlayer[], complex long double *product[], 
	    int, int, int);
    void gauss_jordan(complex long double *surfs[], int, complex long double Cc[]);

    /* 
     *For Tidal Forcing, Surface BCs are:
     *	y3 = 0		Tangential Stress
     *  y4 = 0		Radial Stress
     *  y6 = -(2l+1)/a	Continuity of Potential
     */

     BCs[1] = 0.0;
     BCs[2] = 0.0;
     BCs[3] = (long double)(-(2*L+1)) / layers[max_layers-1].rad;

    for(i=1; i<=6; i++) {
    	Ylayer[i] = (complex long double *)malloc(7*sizeof(complex long double));
    	Yinvlayer[i] = (complex long double *)malloc(7*sizeof(complex long double));
    	product[i] = (complex long double *)malloc(4*sizeof(complex long double));
    	propmat[i] = (complex long double *)malloc(4*sizeof(complex long double));
    }

    for(i=1; i<=3; i++)
    	surfs[i] = (complex long double *)malloc(4*sizeof(complex long double));
    
    /* Form Propagator Matrix */
    /* Start from bottom and work up (new factor mult. from left) */
    /* Set initial Prop. Mat. as Ic Interface Matrix */
    for(i=1;i<=6;i++) {
    	for(j=1;j<=3;j++) {
    		propmat[i][j] = Ic[i][j];
    	}
    }

    for (layer = 1; layer <= (max_layers-1); layer++) {

//	(void)fprintf(stdout,"Layer = %d\n",layer);
    	/* Get Y, Y' for layer */
        for(i=1;i<=6;i++)
        	for(j=1;j<=6;j++) {
        		Ylayer[i][j] = Y[i][j][layer];
        		Yinvlayer[i][j] = Yinv[i][j][layer];
        	}
	
        /* product = Y'(layer) * propmat */
        mat_multi_mat(Yinvlayer,propmat,product,6,3,6);
        for(i=1;i<=6;i++)
        	for(j=1;j<=3;j++)
        		propmat[i][j] = product[i][j];
	
        /* product = Y(layer) * propmat */
        mat_multi_mat(Ylayer,propmat,product,6,3,6);
        	for(i=1;i<=6;i++)
        		for(j=1;j<=3;j++)
        			propmat[i][j] = product[i][j];

    }

    /*
     * Rows 3, 4, and 6 correspond to equations that have BCs at the surface
     * Rewrite these rows in a new matrix surfs*
     *
     * surfs	* Cc	= BCs
     * 3x3	* 3x1	= 1x1
     */
    for(j=1;j<=3;j++) {
    	surfs[1][j] = propmat[3][j];
    	surfs[2][j] = propmat[4][j];
    	surfs[3][j] = propmat[6][j];
    }

    /* Now solve for Cc */
    /* System is simple enough that Gauss-Jordan should be ok */

    for(j=1;j<=3;j++) 	    /*  gauss_jordan replaces input vector with soln */
    	Cc[j] = BCs[j];	    /*  copy into soln vector, so as not to lose it */

    gauss_jordan(surfs,3,Cc);

    /* Free Memory */
    for (i=1; i<=6; i++) {
    	free((complex long double *)Yinvlayer[i]);
    	free((complex long double *)Ylayer[i]);
    	free((complex long double *)propmat[i]);
    	free((complex long double *)product[i]);
    }
    for (i=1; i<=3; i++)
    	free((complex long double *)surfs[i]);

}

/************************************************************
 * gauss_jordan -- Solves a linear system of equations	    *
 *							                                *
 *		A x = b					                            *
 *							                                *
 *		using Gauss-Jordan Elimination		                *
 *							                                *
 * Parameters						                        *
 *		a_mat	-- Input Matrix	(replaced by its            *
 *			    inverse on output)		                    *
 *		b_vec	-- RHS vector (replaced by solution         *
 *			    vector on output		                    *
 *		n	-- number of equations		                    *
 ************************************************************/

void gauss_jordan(complex long double *a_mat[], int n,
		complex long double b_vec[])
{

	int i,j,k;						/* Row and Column Indices */
	int i_max;						/* Row Index of max value in column */

	long double a_max;				/* Max value in column */
	complex long double temp;		/* Variable for swapping */
	complex long double norm;		/* Variable for row reduction */
	complex long double *a_orig[4];	/* Copy of original A matrix */
	complex long double *bsol;		/* b from solution */

	/* Copy a_mat into a new matrix. The original a_mat gets succesively lost
	 * and replaced by its inverse.
	 */

    for(i=1;i<=n;i++) {
        a_orig[i] = (complex long double *)malloc((n+1)*sizeof(complex long double));
        for(j=1;j<=n;j++)
        	a_orig[i][j] = a_mat[i][j];
    	}
    bsol = (complex long double *)malloc(sizeof(complex long double));

/*    fprintf(stderr,"A_orig \t \t \t b:\n");
        for(i=1;i<=n;i++) {
        	for(j=1;j<=n;j++) {
        	    fprintf(stderr,"%20.16Lg + %20.16Lg i\t",creall(a_orig[i][j]),
        	    								cimagl(a_orig[i][j]));
        	    }
        	fprintf(stderr,"| \t %20.16Lg \n",creall(b_vec[i]));
        }*/

	for(k=1;k<=n;k++){

		/* Find the pivot for column k */
		i_max = k;
		a_max = 0.0;
		for(i=k;i<=n;i++) {
			//printf("i = %d, k = %d, %e\n", i, k, creal(a_mat[i][k])); // TODO delete
			if(creal(a_mat[i][k]) >= a_max) {
				i_max = i;
				//printf("Re(a_mat[%d][%d]) = %e\n", i, k, creal(a_mat[i][k])); // TODO delete
				//printf("abs(a_mat[%d][%d]) = %e\n", i, k, cabsl(a_mat[i][k])); // TODO delete
				a_max = a_mat[i][k];
				//printf("a_max set to %e\n", a_max); // TODO delete
			}
		}
		//printf("a_max = %Lf\n", a_max); // TODO delete
		if(a_max == 0.0) {
			fprintf(stderr,"gauss_jordan: Matrix is singular!\n");
		}

		/* Swap rows */
		for (j=1;j<=n;j++)
    		SWAP(a_mat[k][j],a_mat[i_max][j])
    	SWAP(b_vec[k],b_vec[i_max])

    	/* For all rows below pivot */
    	for (i=k+1; i<=n; i++) {
    		/* For all remaining elements in current row */
    		for (j=k+1; j<=n; j++) {
    			a_mat[i][j] -= a_mat[k][j] * (a_mat[i][k]/a_mat[k][k]);
    		}
    		b_vec[i] -= b_vec[k] * (a_mat[i][k]/a_mat[k][k]);
    		/* Fill lower triangular matrix with zeroes */
    		a_mat[i][k] = 0.0;
    	}
	}

	/* Now in Upper Triangular (echelon form)
     * Reduce rows from bottom up
     */
    for (k=n; k>=1; k--) {

    	/* Normalize kth row */
    	norm = a_mat[k][k];
    	for (j=k; j<=n; j++) {
    		a_mat[k][j] /= norm;
    	}
    	b_vec[k] /= norm;

    	/* Reduce overlying rows to get zero in kth column */
    	for (i=(k-1); i>=1; i--) {
    		norm = a_mat[i][k];
    		for (j=k; j<=n; j++) {
    		    a_mat[i][j] -= a_mat[k][j] * norm;
    		}
    		b_vec[i] -= b_vec[k] * norm;
    	}

    }

/*    fprintf(stderr,"A reduced row echelon \t \t \t | \t x:\n");
    for(i=1;i<=n;i++) {
       	for(j=1;j<=n;j++) {
       	    fprintf(stderr,"%7.5Lg + %7.5Lg i\t",creall(a_mat[i][j]),
       	    								cimagl(a_mat[i][j]));
       	    }
       	fprintf(stderr,"| \t %7.5Lg + %7.5Lg i\n",creall(b_vec[i]),
       	        	                       	    		cimagl(b_vec[i]));
    }*/

/*    // Check solution
    for (i=1; i<=n; i++)
   		bsol[i] = 0.0;

    // b = A * x  (b_vec now contains x)
    for (i=1; i<=n; i++)
     	for(j=1; j<=n; j++) {
        	bsol[i] += a_orig[i][j] * b_vec[j];
        }

    fprintf(stderr,"Ax:\n");
        for(i=1;i<=n;i++) {
           	fprintf(stderr,"%Lg + %Lg i\n",creall(bsol[i]),cimagl(bsol[i]));
        }*/

    for(i=1;i<=n;i++) {
    	free((complex long double *)a_orig[i]);
    }
	free(bsol); // DEBUG

}


/************************************************************
 * tidal_potential -- Determine the tidal potential at the  *
 *		surface, varying in lat and long.	                *
 *							                                *
 * Parameters						                        *
 *		t	    	-- time		                            *
 *		potential   -- Tidal Potential		                *
 *		dpot	    --  derivatives of pot	                *
 *				0 = d/dtheta		                        *
 *				1 = d/dphi		                            *
 *				2 = d^2/dtheta^2       		                *
 *				3 = d^2/dphi^2		                        *
 *				4 = d^2/(dtheta dphi)	                    *
 ************************************************************/

void tidal_potential(long double t, long double **potential, 
	long double **dpot[]) {

    int i,j;	/* counters for lat. and long. */
    int n;	/* counter for time */
    int num_time;   /* Number of time intervals */
    long double r,ecc,omega;	/* temp variables */
    long double th,f;		
    char filename[20];	/* Name of output file */
    static int here = 0;    /* Number of times in this function */

    FILE *out_file;

    long double legendre(int,int,long double);

    here++;

    r = grids.radius[grids.num_r-1];
    ecc = planet.ecc;
    omega = planet.freq;

    /* Start at element 1 */
    for(i=0;i<=grids.num_theta;i++) {
    	th = grids.etheta[i];

    	for(j=0;j<=grids.num_phi;j++) {
    		f = grids.ephi[j];

    		potential[i][j] = r*r * omega*omega * ecc
    		    	* (-1.5 * legendre(2,0,cos(th)) * cos(omega*t)
    		    			+ 0.25 * legendre(2,2,cos(th)) * (3.0*cos(2.0*f)
    						* cos(omega*t) + 4.0*sin(2.0*f) * sin(omega*t) ) );

    		dpot[0][i][j] = 3.0 * r*r * omega*omega * ecc * sin(th)*cos(th)
				* cos(f) * (3.0*cos(f) * cos(omega*t)
						+ 4.0*sin(f) * sin(omega*t) );

    		dpot[1][i][j] = 0.75 * r*r * omega*omega * ecc * sin(th)*sin(th)
				* (-7.0*sin(2.0*f-omega*t) + sin(2.0*f+omega*t));

    		dpot[2][i][j] = -1.5 * r*r * omega*omega * ecc * cos(2.0*th)
    			* cos(f) * (-7.0*cos(f-omega*t) + cos(f+omega*t));

    		dpot[3][i][j] = 1.5 * r*r * omega*omega * ecc * sin(th)*sin(th)
				* (-7.0*cos(2.0*f-omega*t) + cos(2.0*f+omega*t));

    		dpot[4][i][j] = -3.0 * r*r * omega*omega * ecc * cos(th) * sin(th)
				* (3.0 * sin(2.0*f) * cos(omega*t)
						- 4.0*cos(2.0*f) * sin(omega*t) );

    	}
    }
}

/************************************************************
 * legendre -- Computes the associated Legendre polynomial  *
 *		P_l^m(x).				                            *
 *							                                *
 * Parameters						                        *
 *		l   --	degree, l >= 0			                    *
 *		m   --	order,	0 <= m <= l		                    *
 *		x   --	cos(theta), -1 <= x <= 1	                *
 *															*
 * Returns													*
 *		plm --  P_l^m(x)									*
 ************************************************************/

long double legendre(int l, int m, long double x)
{
    long double plm,pmm,pmplus1m,plminus1m,plminus2m; // Various Leg. Polyns.
    long double dblfact; // Double factorial
    long double temp;
    int i, ll;

    /* Error checking */
    if (l < 0)
    	fprintf(stderr,"Must have l >= 2\n");
    if (m < 0 || m > l)
    	fprintf(stderr,"Must have 0 <= m <= l\n");
    if (fabs(x) > 1.0)
    	fprintf(stderr,"Must have -1 <= cos(theta) <= 1\n");

    /* Use recurrence relation
     *  (l-m+1) P_(l+1)^m = (2l+1) x P_l^m - (l+m)P_(l-1)^m
     *  (l-m) P_l^m = (2l-1) x P_(l-1)^m - (l-1+m)P_(l-2)^m
     */

    /* Case where l=m has closed form solution:
     * P_l^l = (-1)^l (2l-1)!! (sqrt(1-x^2))^l
     * P_m^m = (-1)^m (2m-1)!! (sqrt(1-x^2))^m
     */
    temp = sqrt(1.0-x*x);	/* This factor shows up a lot */
    pmm = 1.0;					/* P_0^0 = 1

    /* For l=m > 0:
     * Compute double factorial
     */
    dblfact=1.0;
    for (i=1; i<=m; i++)
    	dblfact *= (2.0*(float)i - 1.0);
    /* Account for factors of -1 and temp */
    pmm = dblfact * pow(-temp,m); /* m const in recur rel., index by m */

    /* If l does equal m, we're done. */
    if (l == m)
    	plm = pmm;
    /* Otherwise, use recurrence relation above: */
    else if (l > m) {
    	/* Calculate P_(l+1)^l (or P_(m+1)^m). Note that recurrence relation
    	 * will have a P_(l-1)^(l) term in it, which is identically zero since
    	 * l > m.
    	 */
    	pmplus1m = x*(2*m+1)*pmm;
    	/* If l equal m+1, we're done. */
    	if (l == (m+1))
    		plm = pmplus1m;
    		/* Otherwise, use recurrence relation with both terms: */
    	else {
    		plminus2m = pmm;	/* Set P_(l-2)^m */
    		plminus1m = pmplus1m;		/* Set P_(l-1)^m */
    		for (i=(m+2); i<=l; i++) {
    			plm = (x*(2*i-1)*plminus1m - (i-1+m)*plminus2m) / (i-m);
    			plminus2m = plminus1m;
    			plminus1m = plm;
    		}
    	}
    }

return plm;
}


/************************************************************
 * stress_strain -- Computes the stress and strain at each  *
 *		point within the body resulting from tidal          *
 *		forcing.				                            *
 *							                                *
 * Parameters						                        *
 *		y_i	--  radial functions		                    *
 *		pot	--  tidal potential		                        *
 *		dpot	--  derivatives of pot		                *
 *				0 = d/dtheta		                        *
 *				1 = d/dphi		                            *
 *				2 = d^2/dtheta^2	                        *
 *				3 = d^2/dphi^2		                        *
 *				4 = d^2/(dtheta dphi)	                    *
 *		strain	--  strain tensor		                    *
 *		stress	--  stress tensor		                    *
 *				1 = theta theta		                        *
 *				2 = phi phi		                          	*
 *				3 = r r			                            *
 *				4 = theta phi		                       	*
 *				5 = phi r		                            *
 *				6 = r theta		                            *
 ************************************************************/

void stress_strain(int max_layers, complex long double *y_i[], 
	long double **pot, long double **dpot[], complex long double *strain[],
	complex long double *stress[])
{
    int i,j,k;				/* spatial counters */
    int p;			        /* component counter */
    int surfnode, node;			/* element number */
    int layer;				/* the layer our point is in */
    long double th,f,r;			/* temp vars */
    complex long double lambda, mu;	/* Lame Constants */
    complex long double hydrostatic;	/* Hydrostatic stress */

    FILE *out_file;

    int find_layer(long double, int);

    //out_file = fopen("sstest.dat","w");

    /* At each element */

    for (i=1; i<=(grids.num_theta-1); i++) {
    	th = grids.etheta[i];
    	for (j=0; j<=(grids.num_phi-0); j++) {
    		f = grids.ephi[j];
    		surfnode = i*(grids.num_phi+1) + j;
    		for (k=0; k<=grids.num_r; k++) {
    			r = grids.eradius[k];
    			node = i*(grids.num_phi+1)*(grids.num_r+1) + j*(grids.num_r+1) + k;
    			/* Find which layer we're in and assign mat. properties*/
    			layer = find_layer(r,max_layers);

    			lambda = layers[layer].lambda;
    			mu = layers[layer].rigidity;
    			//(void)fprintf(out_file,"%d %d %d %Lf %Lf %Lf %d %d %.2Le %.2Le %.2Le %.2Le\n",i,j,k,th,f,r,node,layer,creall(lambda),cimagl(lambda),creall(mu),cimagl(mu));

    			/*
    			 * mu = 0 at CMB.  Will cause strain to be undefined
    			 * use Lame params from layer 1 here.
    			 */

    			if(layer==0) {
    				lambda = layers[1].lambda;
    				mu = layers[1].rigidity;
    			}

    			/* Get each component of the strain tensor */

    			strain[1][node] = 2.0/r * pot[i][j] * y_i[1][k]
    			                + 2.0/r * dpot[2][i][j] * y_i[2][k];

    			strain[2][node] = 2.0/r * pot[i][j] * y_i[1][k]
    			                + 2.0/(r*sin(th)) * (1.0/sin(th)
    			                		* dpot[3][i][j]
    			                + cos(th)*dpot[0][i][j]) * y_i[2][k];

    			strain[3][node] = (-4.0*lambda/(lambda+2.0*mu))
				    					* 1.0/r * pot[i][j] * y_i[1][k]
				    			+ ((long double)(2*L*(L+1))
				    					*lambda/(lambda+2.0*mu))
				    					* 1.0/r * pot[i][j] * y_i[2][k]
				    			+ (2.0/(lambda+2.0*mu))
				    					* pot[i][j] * y_i[3][k];

    			strain[4][node] = 2.0/(r*sin(th)) * (dpot[4][i][j]
    			                - 1.0/tan(th) * dpot[1][i][j]) * y_i[2][k];

    			strain[5][node] = 1.0/(mu*sin(th)) * dpot[1][i][j] * y_i[4][k];

    			strain[6][node] = 1.0/mu * dpot[0][i][j] * y_i[4][k];
		
//		(void)fprintf(out_file,"%.2Le %.2Le %.2Le %.2Le %.2Le %.2Le\n",creall(strain[1][node]),creall(strain[2][node]),creall(strain[3][node]),creall(strain[4][node]),creall(strain[5][node]),creall(strain[6][node]));
//		(void)fprintf(out_file,"%.2Le %.2Le %.2Le %.2Le %.2Le %.2Le\n",cimagl(strain[1][node]),cimagl(strain[2][node]),cimagl(strain[3][node]),cimagl(strain[4][node]),cimagl(strain[5][node]),cimagl(strain[6][node]));

    			for(p=1;p<=6;p++)
    				strain[p][node] *= 0.5;

    			/* Convert Stresses into Strains */
    			/* Hydrostatic stress */
    			hydrostatic = lambda*(strain[1][node] + strain[2][node]
    			                                      + strain[3][node]);

    			for(p=1;p<=6;p++)
    				stress[p][node] = 2.0*mu*strain[p][node];

    			/* Add hydro. stress to diagonal terms */
    			for(p=1;p<=3;p++)
    				stress[p][node] += hydrostatic;

    			/* Cannot solve for all components at pole */
    			/* Use strain at adjacent theta */
    			if(i==1) {
    				for(p=1; p<=6; p++) {
    					strain[p][node - (grids.num_phi+1)*(grids.num_r+1)]
    					          = strain[p][node];
    					stress[p][node - (grids.num_phi+1)*(grids.num_r+1)]
    					          = stress[p][node];
    				}
    			}
    			if(i==grids.num_theta-1) {
    				for(p=1; p<=6; p++) {
    					strain[p][node + (grids.num_phi+1)*(grids.num_r+1)]
    					          = strain[p][node];
    					stress[p][node + (grids.num_phi+1)*(grids.num_r+1)]
    					          = stress[p][node];
    				}
    			}

    		}
    	}
    }

    // (void)fclose(out_file);

}

/************************************************************
 * find_layer -- Determines which layer a radial point is   *
 *		in.				 	                                *
 *							                                *
 * Parameters						                        *
 *		r   --	radial position of point	                *
 *							                                *
 * Returns						                            *
 *		layer	--  layer the point is in	                *
 ************************************************************/

int find_layer(long double r, int max_layers) {

    int k;

    for(k=0; k<=(max_layers-1); k++) {
    	if (r < layers[k].rad) {
    		return(k);
    		break;
    	}
    }

    /* If we're at the surface, treat that as layer below */
    	return(max_layers-1);
}


/************************************************************
 * heating -- Computes the tidal heating in each element    *
 *							                                *
 * Parameters						                        *
 *		strain   --	strain			                        *
 *		stress	 --	stress			                        *
 *		tidal_heat --	volumetric heating rate	            *
 *							                                *
 ************************************************************/

void heating(complex long double *strain[], complex long double *stress[],
	long double *tidal_heat) {

    int i,j,k;			/* positional counters */
    int p;			/* component counter */
    int node, surfnode;		/* Nodal indices */
    static int been_here = 0;
    char filename[80];

    FILE *out_file;

/*  sprintf(filename,"heattest.%d.dat",been_here);
    out_file = fopen(filename,"w");
*/
    for (i=0; i<=grids.num_theta; i++) {
    	for (j=0; j<=grids.num_phi; j++) {
    		surfnode = i*(grids.num_phi+1) + j;
    		for (k=1; k<=grids.num_r; k++) {
    			node = i*(grids.num_phi+1)*(grids.num_r+1) + j*(grids.num_r+1) + k;

    			/*
    			 * If this is the first time through, initialize
    			 * the heating to 0 everywhere, before adding in
    			 * the contributions from the timestep.
    			 */
    			if(been_here == 0)
    				tidal_heat[node] = 0.0;

    			for (p=1;p<=3;p++) {
    				tidal_heat[node] += planet.freq / NUM_T *
    						(cimagl(stress[p][node])*creall(strain[p][node])
    								- creall(stress[p][node])*cimagl(strain[p][node]));
    			}

    			/* Add 2x off-diagonal components */
    			for (p=4;p<=6;p++) {
    				tidal_heat[node] += planet.freq / NUM_T * 2.0 *
    						(cimagl(stress[p][node])*creall(strain[p][node])
    								- creall(stress[p][node])*cimagl(strain[p][node]));
    			}

    		}
    	}
    }

    been_here++;
//    (void)fclose(out_file);

}

/************************************************************
 * output_data -- Writes relevant data to files.            *
 *		May modify input file to have flags	                *
 *		specifying which fields are written.	            *
 *							                                *
 *		Also does some averaging needed for some            *
 *		data sets.				                            *
 *							                                *
 * Parameters						                        *
 *		max_layers  --	Total layers in planet	            *
 *		y_i	    	--	Radial functions	                *
 *		potential   --	Tidal Potential		                *
 *		dpot	    --	1st and 2nd derivatives of          *
 *						tidal potential		                *
 *		strain	    --	strain tensor		                *
 *		stress	    --	stress tensor		                *
 *		tidal_heat  --	tidal heating rate	                *
 *		output_flags	--  which files to output           *
 *		directory   --	output directory	                *
 *		tc_iter	    --  which iteration this is             *
 ************************************************************/

void output_data(int max_layers, complex long double *y_i[], 
	long double **potential, long double **dpot[], 
	complex long double *strain[], complex long double *stress[], 
	long double *tidal_heat, char output_flags, char directory[], 
	int tc_iter) {

    int i,j,k;			    			/* positional counters */
    int node, surfnode, axinode;    	/* nodal index */
    int maxnodes, maxsurf, maxaxi;  	/* Total number of nodes */
    int p;			    				/* component counter */
    int layer;
    char filename0[80];		    		/* filename */
    char filename[80], filename1[80];	/* filename */
    long double *heating_axi;	    	/* phi-avged heating, for axisymmetric */
    long double *heating_radial;		/* radially-avged heating */
    long double dr, dth, dvol, dvol2;	/* Volume elements */
    long double heating_ave;	    	/* Global heating total */
    long double have_sil, have_ice; 	/* Global heating layers */
    long double hsurf_ave,hbot_ave; 	/* Average heating at surfaces */
    long double vol_ice, vol_sil;   	/* Volumes */
    long double vol_surf,dsurf;			/* Surface */
    long double Hmu,hHmu;	    		/* Sensitivity Parameter */
    long double r,K;		    		/* temp params */
    complex long double mu,dy1dr;
    complex long double Loveh2, Lovek2;	/* Love Number */
    long double phaseh2, phasek2;		/* Phase Lags */

    FILE *out_file0;					/* output file-stream */
    FILE *out_file;						/* output file-stream */
    FILE *out_file2;					/* output file-stream */
    FILE *out_file3;					/* output file-stream */

    int find_layer(long double, int);

    maxnodes = (grids.num_r+1)*(grids.num_theta+1)*(grids.num_phi+1);
    maxsurf = (grids.num_theta+1)*(grids.num_phi);
    maxaxi = (grids.num_theta+1)*(grids.num_r+1);

    heating_axi = (long double *)malloc(maxaxi*sizeof(long double));
    heating_radial = (long double *)malloc((grids.num_r+1)*sizeof(long double));

    heating_ave = 0.0;
    hsurf_ave = 0.0;
    hbot_ave = 0.0;
    have_sil = 0.0;
    have_ice = 0.0;
    vol_sil = 0.0;
    vol_ice = 0.0;
    vol_surf = 0.0;

    /* out_file0 for tests and debugging.  Don't use for real data */
    /*sprintf(filename0,"%s/test_output.dat",directory);
    out_file0 = fopen(filename0,"a");*/
    /*fprintf(stderr,"%s\n",filename0);*/

/*    (void)fprintf(stderr,"\tOutputs %x\n",output_flags);*/
/*    directory[strlen(directory)-1] = '\0';*/
/*    (void)fprintf(stderr,"\tOutputs %x\n",output_flags);*/
    //(void)fprintf(stderr,"directory %s\n",directory);

    /* Get axisymmetric, radial, and average heating */
	for (k=0; k<=(grids.num_r); k++) {	/* Don't use surface node */

	  if(k==0)
          dr = grids.radius[0];   
      else if (k==grids.num_r) {
          dr = 0.0;
          dvol2 = 0.0;
      }
	  else 
		  dr = (grids.radius[k] - grids.radius[k-1]);
	
		/* Initialize */
		heating_radial[k] = 0;

    for (i=0; i<=(grids.num_theta-0); i++) {

/*	if(i==1)
	    dth = (theta[i+1] - theta[i]) / 2.0;
	else if(i==grids.num_theta)
	    dth = (theta[i] - theta[i-1]) / 2.0;
	else 
	    dth = (theta[i+1] - theta[i-1]) / 2.0;
*/
/*	dth = (grids.theta[i] - grids.theta[i-1]);*/
	dth = grids.dth;
	if (i == 0 || i == grids.num_theta)
	    dth = 0.0;

    dsurf = 2.0*PI*sin(grids.etheta[i])*dth;				/* Nondim surface area element */
		dvol = dsurf*grids.eradius[k]*grids.eradius[k]*dr;	/* Volume of shell */
        if (k<grids.num_r)
          dvol2 = dsurf/3.0 * (pow(grids.radius[k],3) - pow(grids.radius[k-1],3))
;
	   /* if(i==4)
	    (void)fprintf(stderr,"dvol %d %d : %Lg %Lg %Lg %Lg = %Lg\n",
		    i,k,dth,dr,grids.etheta[i],grids.eradius[k],dvol);
*/
	    axinode = i*(grids.num_r+1) + k;
	    heating_axi[axinode] = 0.0;	    /* Initialize */
	    
	    /* Add up all heating at each longitude */
	    /* Don't double-count 0 and 2 pi */
	    for (j=1; j<=(grids.num_phi-1); j++) {
				node = i*(grids.num_phi+1)*(grids.num_r+1) + j*(grids.num_r+1) + k;
				heating_axi[axinode] += tidal_heat[node];
	    }
	    heating_axi[axinode] /= grids.num_phi;	/* Normalize, dphi constant */

		heating_radial[k] += heating_axi[axinode]*dsurf; /* Normalize by surf area of element*/
	    heating_ave += heating_axi[axinode] * dvol;

	    /* Get avg in ice layer */
	    if (k >= 3) {
				heating_axi_ice[i][k-2] = heating_axi[axinode];
				have_ice += heating_axi[axinode] * dvol;
				vol_ice += dvol;
       /*         if (i==grids.num_theta) {
                 fprintf(stdout,"ice %d %Lg %Lg %Lg\n",k,grids.radius[k],
                                        dr,vol_ice); 
                }
                */
	    }

	    /* Get avg in silicate layer */
	    if (grids.radius[k] <= layers[1].rad) {
				have_sil += heating_axi[axinode] * dvol;
				vol_sil += dvol2;
                /*
                if (i==grids.num_theta) {
                    fprintf(stdout,"sil %d %Lg %Lg %Lg\n",k,grids.radius[k],
                                        dr,vol_sil); 
                }
                */
	    }
	    
	    /* Get avg at surface */
	    if (k == grids.num_r) {
				hsurf_ave += heating_axi[axinode] * dsurf;
				vol_surf += dsurf;
	    }
	    
	    /* Get avg at CMB */
	    if (k == 11) {
				hbot_ave += heating_axi[axinode] * dsurf;
	    }
	    
		}
	  /* Have to normalize radial heating by surface area of sphere */
		heating_radial[k] /= (4.0*PI);

	}
    (void)fprintf(stdout,"htotal = %Lg, havg = %Lg, vol = %Lg\n",heating_ave,heating_ave/planet.vol,planet.vol);
    (void)fprintf(stdout,"hsil= %Lg, havg_sil = %Lg, vol_sil = %Lg\n",have_sil,have_sil/vol_sil, vol_sil);
    (void)fprintf(stdout,"hice= %Lg, havg_ice = %Lg, vol_ice = %Lg\n",have_ice,have_ice/vol_ice, vol_ice);
    (void)fprintf(stdout,"hsurf= %Lg, havg_surf = %Lg, vol_surf = %Lg\n",hsurf_ave,hsurf_ave/vol_surf, vol_surf);
    (void)fprintf(stdout,"hbot= %Lg, havg_bot = %Lg r_bot = %Lg\n",hbot_ave,hbot_ave/vol_surf,grids.radius[3]);
    heating_ave /= planet.vol;

	char* dir = directory;
	if (directory[0] == '\"') { dir++; }

    /* Output full 3D heating */
    if(output_flags & (1<<0)) {
    	sprintf(filename,"%s/heating_3D.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	sprintf(filename,"%s/initial_heating.dat",dir);
    	out_file2 = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	sprintf(filename,"%s/a.initial_heating.dat",dir);
    	out_file3 = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	/* Print coords to a.initial_heating file */
    	for (i=0; i<=(grids.num_theta); i++)
    		fprintf(out_file3,"%Lg\n",grids.etheta[i]);
    	for (j=0; j<=(grids.num_phi); j++)
    		fprintf(out_file3,"%Lg\n",grids.ephi[j]);
    	for (k=3; k<=(grids.num_r-0); k++)
    		fprintf(out_file3,"%Lg\n",grids.eradius[k]/grids.radius[max_layers-1]);

    	for (i=0; i<=(grids.num_theta); i++) {
    		for (j=0; j<=grids.num_phi; j++) {
    			for (k=1; k<=(grids.num_r-0); k++) {
    				node = i*(grids.num_phi+1)*(grids.num_r+1)
    						+ j*(grids.num_r+1) + k;
    				fprintf(out_file,"%Lg %Lg %Lg %Lg\n",grids.etheta[i],
    						grids.ephi[j],grids.eradius[k],tidal_heat[node]);
    				if(grids.eradius[k] >= layers[2].rad) {
    					fprintf(out_file2,"%Lg %Lg %Lg %Lg\n",grids.etheta[i],
    							grids.ephi[j],
    							grids.eradius[k]/grids.radius[max_layers-1],
    							tidal_heat[node]);
	    			fprintf(out_file3,"%Lg\n",tidal_heat[node]);
    				}
    			}
    		}
    	}
	
    	(void)fclose(out_file);
    	(void)fclose(out_file2);
    }

    /* Write surface data */
    if(output_flags & (1<<1)) {
    	sprintf(filename,"%s/heating_surf.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	for (i=1; i<(grids.num_theta); i++) {
    		for (j=1; j<(grids.num_phi); j++) {
    			surfnode = i*(grids.num_phi+1) + j;
    			fprintf(out_file,"%Lg %Lg %Lg %d %d\n",grids.etheta[i],
    					grids.ephi[j],
    					tidal_heat[surfnode*(grids.num_r+1) + grids.num_r],
    					surfnode,surfnode*(grids.num_r+1) + grids.num_r);
    		}
    	}

    	(void)fclose(out_file);

    	sprintf(filename,"%s/heating_bot.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	for (i=1; i<(grids.num_theta-1); i++) {
    		for (j=1; j<(grids.num_phi-1); j++) {
    			surfnode = i*(grids.num_phi+1) + j;
    			fprintf(out_file,"%Lg %Lg %Lg\n",grids.etheta[i],grids.ephi[j],
    					tidal_heat[surfnode*(grids.num_r+1) + 11]);
    		}
    	}

    	(void)fclose(out_file);
    }

    /* Write Axisymmetric Heating */
    if(output_flags & (1<<2)) {
    	sprintf(filename,"%s/heating_axi.%d.dat",dir,tc_iter);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	sprintf(filename,"%s/heating_axi_ice.%d.dat",dir,tc_iter);
    	out_file2 = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

		// Create heap structs to store info written into files
		//int h_axi_int_p[3] = { sizeof(grids.etheta) / sizeof(grids.etheta[0]), sizeof(grids.eradius) / sizeof(grids.eradius[0]), sizeof(heating_axi) / sizeof(heating_axi[0])};
		//Tiradata* h_axi = construct_Tiradata("heating_axi", h_axi_int_p);
		//for (i = 0; i < (grids.num_theta + 1) * (grids.num_r + 1); i++) {} // TODO: FILL STRUCT

        for (i=0; i<=(grids.num_theta-0); i++) {
        	for (k=0; k<=(grids.num_r-0); k++) {
        		axinode = i*(grids.num_r+1) + k;
        		fprintf(out_file,"%Lg %Lg %Lg\n",
        				grids.etheta[i],grids.eradius[k],heating_axi[axinode]);

        		/* Print out only data in ice shell */
        		/* if single layer used for silicate core */
        		if (k >= 3) {
        			fprintf(out_file2,"%Lg %Lg %Lg\n",
        					grids.etheta[i],grids.eradius[k],heating_axi[axinode]);
        		}
        	}
        }
        (void)fclose(out_file);
        (void)fclose(out_file2);
    }
    
    /* Write Radial Heating */
    if(output_flags & (1<<3)) {
		sprintf(filename,"%s/heating_radial.%d.dat",dir,tc_iter);
		out_file = fopen(filename,"w");
		fprintf(stdout,"%s\n",filename);

		sprintf(filename,"%s/heating_radial_ice.%d.dat",dir,tc_iter);
		out_file2 = fopen(filename,"w");
		fprintf(stdout,"%s\n",filename);

		for (k=0; k<=(grids.num_r-0); k++) {
		    fprintf(out_file,"%Lg %Lg\n", grids.eradius[k],heating_radial[k]);

		    /* Print out only data in ice shell */
		    if (k >= 3) {
				fprintf(out_file2,"%Lg %Lg\n",grids.eradius[k],heating_radial[k]);
		    }
	    }	

		(void)fclose(out_file);
		(void)fclose(out_file2);
    }
    
    /* Write potential */
    if(output_flags & (1<<4)) {
    	sprintf(filename,"%s/potential.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

        for (i=0; i<=(grids.num_theta-0); i++) {
        	for (j=0; j<=(grids.num_phi-0); j++) {
        		fprintf(out_file,"%Lg %Lg %Lg\n",grids.etheta[i],grids.ephi[j],
        				potential[i][j]);
        	}
        }

        (void)fclose(out_file);
    }

    /* Write Radial functions. */
    if(output_flags & (1<<5)) {
    	fprintf(stderr,"%s\n",dir);
    	sprintf(filename,"%s/y_i_real.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);
    	sprintf(filename,"%s/y_i_imag.dat",dir);
    	out_file2 = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	for(k=0;k<=grids.num_r;k++) {
    		(void)fprintf(out_file,"%Lg\t",grids.eradius[k]);
    		(void)fprintf(out_file2,"%Lg\t",grids.eradius[k]);
    		for(i=1;i<=6;i++) {
    			(void)fprintf(out_file,"%Lg\t",creall(y_i[i][k]));
    			(void)fprintf(out_file2,"%Lg\t",cimagl(y_i[i][k]));
    		}
    		(void)fprintf(out_file,"\n");
    		(void)fprintf(out_file2,"\n");
    	}

    	fclose(out_file);
    	fclose(out_file2);

    	/* Output Love Number */
    	Lovek2 = (-y_i[5][grids.num_r] - 1.0);
    	Loveh2 = y_i[1][grids.num_r] * layers[max_layers-1].g;
    	phaseh2 = atan(-cimagl(Loveh2)/creall(Loveh2)) * 180.0 / PI;
    	phasek2 = atan(-cimagl(Lovek2)/creall(Lovek2)) *180.0 / PI;
	
    	sprintf(filename,"%s/lovenumber.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);
    	(void)fprintf(out_file,"h2 = %Lg + %Lg, amplitude %Lg, phase %Lg\n",
    			creall(Loveh2),cimagl(Loveh2),cabsl(Loveh2),phaseh2);
    	(void)fprintf(out_file,"k2 = %Lg + %Lg, amplitude %Lg, phase %Lg\n",
    			creall(Lovek2),cimagl(Lovek2),cabsl(Lovek2),phasek2);
    	fclose(out_file);
    
    	sprintf(filename,"%s/Imk2.dat",dir);
    	out_file2 = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);
    	(void)fprintf(out_file2,"%Lg\n",cimagl(Lovek2));
    	fclose(out_file2);

    	/* Tobie's Sensitivity Parameters */
    	sprintf(filename,"%s/sensparam.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);
    	for(k=0;k<=grids.num_r;k++) {
    		r = grids.eradius[k];
    		layer = find_layer(r,max_layers);
    		K = layers[layer].bulk;
    		mu = layers[layer].rigidity;
    		dy1dr = conjl(-2.0/r * y_i[1][k]
    		                              + (long double)(L*(L+1))/r * y_i[2][k]);
    		Hmu = ( (4.0/3.0) * (r*r)/pow(cabs(K + (4.0/3.0)*mu),2)
    				* powl(cabsl(y_i[3][k] - cabsl(K + (4.0/3.0)*mu)/r
    				     * (2.0*y_i[1][k] - (long double)(L*(L+1))*y_i[2][k])),2)
    			- (4.0/3.0)*r * creall(dy1dr*(2.0*y_i[1][k]
    			                   - (long double)(L*(L+1))*y_i[2][k]))
    			+ (1.0/3.0)*(powl(cabsl(2.0*y_i[1][k]
    			                   - (long double)(L*(L+1))*y_i[2][k]),2))
    			+ (long double)(L*(L+1))*r*r
    			         * (powl(cabsl(y_i[4][k]),2))/(powl(cabsl(mu),2))
    			+ (long double)(L*(L*L-1)*(L+2))*(powl(cabsl(y_i[2][k]),2)) );
    		hHmu = 2.1 * (powl(planet.freq,5) * powl(grids.r_outer,4) * planet.ecc*planet.ecc)/(r*r) * Hmu * cimagl(mu);
    		(void)fprintf(out_file,"%Lg\t",grids.eradius[k]);
    		(void)fprintf(out_file,"%Lg %Lg %Lg %Lg\n",Hmu,hHmu,creall(mu),cimagl(mu));
    	}
    	fclose(out_file);

    }
  
    /* Write Stress and Strain Tensors */
    if(output_flags & (1<<6)) {
    	sprintf(filename,"%s/stress_real.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);
    	sprintf(filename,"%s/stress_imag.dat",dir);
    	out_file2 = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	for (i=0; i<=(grids.num_theta-0); i++)
    		for (j=0; j<=(grids.num_phi-0); j++)
    			for (k=1; k<=grids.num_r; k++) {
    				node = i*(grids.num_phi+1)*(grids.num_r+1) + j*(grids.num_r+1) + k;
    				fprintf(out_file,"%.4Lf %.4Lf %4Lf ",
    						grids.etheta[i],grids.ephi[j],grids.eradius[k]);
    				fprintf(out_file2,"%.4Lf %.4Lf %4Lf \n",
    						grids.etheta[i],grids.ephi[j],grids.eradius[k]);
    				for(p=1;p<=3;p++) {
    					fprintf(out_file,"%.6Lg ",creall(stress[p][node]));
    					fprintf(out_file2,"%.6Lg ",cimagl(stress[p][node]));
    				}
    				fprintf(out_file,"%d\n",node);
    				fprintf(out_file2,"%d\n",node);
    			}

    	fclose(out_file);
    	fclose(out_file2);

    	sprintf(filename,"%s/strain_real.dat",dir);
    	out_file = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);
    	sprintf(filename,"%s/strain_imag.dat",dir);
    	out_file2 = fopen(filename,"w");
    	fprintf(stdout,"%s\n",filename);

    	for (i=0; i<=(grids.num_theta-0); i++)
    		for (j=0; j<=(grids.num_phi-0); j++)
    			for (k=1; k<=grids.num_r; k++) {
    				node = i*(grids.num_phi+1)*(grids.num_r+1) + j*(grids.num_r+1) + k;
    				fprintf(out_file,"%.4Lf %.4Lf %4Lf ",
    						grids.etheta[i],grids.ephi[j],grids.eradius[k]);
    				fprintf(out_file2,"%.4Lf %.4Lf %4Lf \n",
    						grids.etheta[i],grids.ephi[j],grids.eradius[k]);
    				for(p=1;p<=3;p++) {
    					fprintf(out_file,"%.6Lg ",creall(strain[p][node]));
    					fprintf(out_file2,"%.6Lg ",cimagl(strain[p][node]));
    				}
    				fprintf(out_file,"%d\n",node);
    				fprintf(out_file2,"%d\n",node);
    			}

    	fclose(out_file);
    	fclose(out_file2);
    }

    free((long double *)heating_axi);
/*    (void)fclose(out_file0);*/
}

long double calculate_c_s(char C_or_S, long double alpha) {
	if (C_or_S == 'C') {
		return lgammal(1 + alpha) * cos(PI * alpha / 2.0);
	} else if (C_or_S == 'S') {
		return lgammal(1 + alpha) * sin(PI * alpha / 2.0);
	} else {
		return 0.0;
	}
}
