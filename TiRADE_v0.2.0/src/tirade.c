/********************************************************
 * tirade -- A program to compute the		            *
 *	conductive temperature structure, and		        *
 *	corresponding tidal heating distribution for a	    *
 *	spherically symmetric planet or moon composed      	*
 *	of an arbitrary number of layers.		            *
 *							                            *
 * Usage:						                        *
 *	tirade.x input_file				                    *
 *							                            *
 *	The user provides an input file containing a	    *
 *	the mechanical properties for each layer in the	    *
 *	body, its orbital and global properties, and	    *
 *	information about which output is needed.  	        *
 *							                            *
 *	The program computes the tidal heating		        *
 *	distribution using the tidal_heating code, then     * 
 *	computes the conductive temperature structure	    *
 *	based on that heating, adjusting the shell        	*
 *	thickness in order to get an equilibrium	        *
 *	solution, using the conduction code.		        *
 *							                            *
 *	The viscosity of each layer is updated based on	    *
 *	the temperature, and the tidal heating is	        *
 *	recomputed based on the new viscosity and shell	    *
 *	thickness.  The temperature is computed from	    *
 *	the new heating.  The program iterates until	    *
 *	it 1) converges to a solution, 2) the shell	        *
 *	melts, or 3) the ocean freezes.			            *
 *							                            *
 * Purpose:						                        *
 *	To compute the distribution of tidal heating in a   *
 *	multilayered spherically symmetric planetary body.	*
 *							                            *
 * Author:  		James H. Roberts				    *
 *							                            *
 * Last Modified:  31 January 2013					    *
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

long double **heating_axi_ice;
long double **T;

struct grid_info grids;
struct globals planet;
struct boundaries bc;
struct matprop *layers;
struct andrade_params andrade;

int main(int argc, char *argv[])
{

    int i,j,k;				/* counters */
    int max_layers;			/* total number of material layers */
    int tc_iter;			/* Tidal/Conduction counter */
    long double misfit, prev_ri;	/* Misfit */
    char output_flags;			/* Controls which files are output */
    char cond_flag;			/* Controls if conduction is run */
    char true_fluid_flag, rheology_flag;
    unsigned char path_flags_found = OFFBYTE; /* Indicates which variables are time-dependent */
    long double* td_ptr[5]; // array of 5 long double pointers
    char directory[80];			/* output directory */
    char filename[80];			/* Test File name */

    FILE *out_file1;
    FILE *out_file2;

    void read_inputs(char *, int *, char *, char *, char*, char*, char*, long double**, char directory[]);
    void set_grids(int, char*);
    void print_stuff(char directory[]);
    void tidal_heating(int, char, char directory[], int, char*);
    void conduction(char directory[], int);
    void compute_visc(int, long double **, char directory[], int);
    long double calculate_c_s(char, long double);

    if (argc < 2) {
	(void) fprintf(stderr, "Usage is: tirade.x <input_file>\n");
	exit(8);
    }
    
    /* Read in data from input file */

    read_inputs(argv[1], &max_layers, &output_flags, &cond_flag, &true_fluid_flag, &rheology_flag, &path_flags_found, td_ptr, directory);
    set_grids(max_layers, &rheology_flag);

    /* Allocate memory for global heating and temperature arrays */
    heating_axi_ice = (long double **)malloc((grids.num_theta+1)
	    *sizeof(long double));
    heating_axi_ice[0] = (long double *)malloc((grids.num_theta+1)
	    *(grids.num_r+1)*sizeof(long double));
    for(i=1;i<=grids.num_theta;i++)
	heating_axi_ice[i] = heating_axi_ice[i-1] + (grids.num_r+1);

    T = (long double **)malloc((grids.num_theta+1)*sizeof(long double));
    T[0] = (long double *)malloc((grids.num_theta+1)*(grids.num_r+1)
	    *sizeof(long double));
    for(i=1;i<=grids.num_theta;i++)
	T[i] = T[i-1] + (grids.num_r+1);

    if (rheology_flag & ANDRADE_RHEOLOGY) {
        andrade.C = calculate_c_s('C', andrade.alpha);
        andrade.S = calculate_c_s('S', andrade.alpha);
    }

    if (cond_flag & (1<<0)) { 
        tc_iter = 0;
        misfit = 1.0;
        while (tc_iter <=1 || (tc_iter <=20 && fabs(misfit) > 0.01) ) {
            prev_ri = grids.r_inner;
            /* Solve for tidal heating */
            tidal_heating(max_layers, output_flags, directory, tc_iter, &true_fluid_flag);

            /* Solve for conductive temperature profile */
            conduction(directory, tc_iter);
            misfit = (grids.r_inner - prev_ri) / prev_ri;

            /* Determine viscosity */
            compute_visc(max_layers, T, directory, tc_iter);

            fprintf(stderr,"*** TC %d = %Le, ri = %Lf ***\n",tc_iter,misfit,grids.r_inner);
            tc_iter++;
        }
    }

    else 
	tidal_heating(max_layers, output_flags, directory, 0, &true_fluid_flag);

    (void)free((struct matprop *)layers);

}

/************************************************************
 * read_inputs -- reads the material and global properties  *
 *		from a user-supplied input file.	                *
 *							                                *
 * Parameters						                        *
 *	input_filename	-- name of the input file	            *
 *	max_layers	-- the total number of layers	            *
 *	output_flags	-- bitwise mask for which files	        *
 *			    should be written to output	                *
 *	cond_flag	-- flag for conduction		                *
 *	directory	-- the output directory		                *
 ************************************************************/

void read_inputs(char *input_filename, int *max_layers_ptr, 
	char *output_flags_ptr, char *cond_flag, char* true_fluid_flag, char* rheology_flag, char* path_flags_found, long double** td_ptr, char directory[]) {

    FILE *in_file;  /* Input file */
    char line[120];  /* line from input file */
    char dir[80];
    int layer;	    /* the current layer */
    int itemp1,itemp2,itemp3;
    int i,j,k;	    /* counters */
    int flag;	    /* the current output flag */
    static int max_layers; /* the total number of layers */
    /* temporary variables from read in */
    long double temp1,temp2,temp3,temp4,temp5;
    long double dr;	    /* radial spacing of nodes */

    in_file = fopen(input_filename,"r");
    fprintf(stdout,"\tReading file %s\n",input_filename);

    /* Read in grid spacing */
    (void)fgets(line,sizeof(line),in_file);
    (void)sscanf(line,"%d %d %d",&itemp1,&itemp2,&itemp3);
    grids.num_theta = itemp1;
    grids.num_phi = itemp2;
    grids.num_r = itemp3;

    /* Read in global data */
    (void)fgets(line,sizeof(line),in_file);
    (void)sscanf(line,"%Le %Le %Lf",&temp1,&temp2,&temp3);
    planet.visc0 = temp1;
    planet.E_a = temp2;
    planet.k = temp3;

    (void)fgets(line,sizeof(line),in_file);
    (void)sscanf(line,"%Le %Le %Le",&temp1,&temp2,&temp3);
    planet.freq = temp1;
    planet.ecc = temp2;
    planet.ob = temp3;

    /* Read in boundary conditions */
    (void)fgets(line,sizeof(line),in_file);
    (void)sscanf(line,"%Lf %Lf %Lf",&temp1,&temp2,&temp3);
    bc.bottemp = temp1;
    bc.surftemp = temp2;
    bc.flux = temp3;
//    (void)fprintf(stderr,"%Lf %Lf %Lf\n",bc.surftemp,bc.bottemp,bc.flux);

    /* Read in output information */
    (void)fgets(line,sizeof(line),in_file);
    //(void)sscanf(line,"%s",directory);
    strcpy(directory, line); // swapped this line out for the one above it so we can get directory names with spaces
    directory = strtok(directory, "\""); // skip to first quote mark wrapping directory
    directory = strtok(directory, "#");
    if (!directory) (void)sscanf(line,"%s",directory); // if quotes not used then just use the line
    (void)fprintf(stderr,"Output directory -> '%s'\n",directory);

    *output_flags_ptr = 0;
    /*
     * Read in flags for output files
     * A 1 in the input file means this file will be written
     * _i_  _field_
     * 0    3D	heating
     * 1    surface heating
     * 2    2D axisymmetric heating
		 * 3		1D radial heating
		 * 4    potential
     * 5    y_i
     * 6    Stress and strain
     */
    char* output_types[30]={"3D heating", "Surface heating", "2D axisymmetric heating", "1D radial heating", "Potential", "Radial functions", "Stress & strain"};
    for(i=0;i<=6;i++) {
        (void)fgets(line,sizeof(line),in_file);
        (void)sscanf(line,"%d",&flag);
        if (flag) {
            *output_flags_ptr |= (1<<i);
            if (*output_flags_ptr && !(*output_flags_ptr & (*output_flags_ptr - 1))) { (void)fprintf(stderr,"\n\t Outputs -> %s\n", output_types[i]); }
            else { (void)fprintf(stderr,"\t         -> %s\n", output_types[i]); }
        }
    }
    if (*output_flags_ptr) { fprintf(stderr,"\n"); }
    //(void)fprintf(stderr,"\tOutputs %x\n",*output_flags_ptr);

    /* Read in conduction flag */
    *cond_flag = 0;
    (void)fgets(line,sizeof(line),in_file);
    (void)sscanf(line,"%d",&flag);
    if (flag) {
	    *cond_flag |= (1<<0);
    }

    /* Read in true fluid flag */
    *true_fluid_flag = 0;
    (void)fgets(line,sizeof(line),in_file);
    (void)sscanf(line,"%d",&flag);
    if (flag) {
	    *true_fluid_flag |= (1<<0);
        // Get cutoff values if a true fluid evaluation is requested
        (void)sscanf(line,"%d %Le %Le",&flag, &planet.true_fluid_cutoff_rigidity, &planet.true_fluid_cutoff_visc);
    }

    /* Read in Maxwell or Andrade rheology flag... can be extended for more rheologies */
    *rheology_flag = 0;
    (void)fgets(line,sizeof(line),in_file);
    (void)sscanf(line,"%d",&flag);
    if (flag==ANDRADE_RHEOLOGY) {
	    *rheology_flag |= ANDRADE_RHEOLOGY;
        // Get alpha and beta
        (void)sscanf(line,"%d %Lf %Le",&flag, &andrade.alpha, &andrade.beta);
        if (andrade.alpha >= ANDRADE_ALPHA_TOO_LARGE) {
            // Approach Maxwell solution (avoid computation errors this way)
            andrade.alpha = 1;
            andrade.beta = 0;
        }
    }

    /* Read in layer information */
    (void)fgets(line,sizeof(line),in_file);
    (void)sscanf(line,"%d",max_layers_ptr);

    max_layers = *max_layers_ptr;
    fprintf(stdout,"\t No. layers -> %d\n\n",max_layers);
    layers = (struct matprop *)malloc((max_layers+2)*sizeof(struct matprop));

    int numFlags = 0;
    char flags[5][4] = { "rad", "den", "she", "bul", "vis" }; // flags to find
    (void)fgets(line,sizeof(line),in_file); // collect flag OR 0th-layer line
    while (flagSearch(line) != OFFBYTE) { // while flags are present in line
        // collect number of time data points they want (will collect numpoints + 1 for label) **
        // collect long double pointer from readCSVData(numRows, numLayers, tdFile) **
        char flaggedLine[80], tdPathname[80];
        char flippedFlag = flagSearch(line); // create mask
        for (i = 0; i < 8; i++) { flippedFlag ^= 1UL << i; } // toggle bits
        *path_flags_found |= flagSearch(line); // indicate which flags have been found

        int l, flagIndex;
        for (l = 0; l < 5; l++) {
            if (strstr(line, flags[l]) != NULL) {
                strcpy(flaggedLine, (strstr(line, flags[l]) + 4));
                flagIndex = l;
                break;
            }
        }
        // say our flaggedLine is "pathname,x" -> parse into tdPathname string & atoi numRows
        char* flagSplit_ptr = strtok(flaggedLine, ","); // parse string
        strcpy(tdPathname, flagSplit_ptr); // copy first parsed part to pathname string
        flagSplit_ptr = strtok(NULL, ","); // go to next string (number of td-data in string format)
        int numRows = atoi(flagSplit_ptr); // collect number of rows (excluding column headers) that's read
            // NOTE: warn user of writing more information in the same line
        // throw error if bad pathname or if atoi numRows is invalid

        // Open requested file, read data with relevant parameters
        // change flagging to give
        printf("\nTime-dependent %s data detected at path: %s, with number of data points requested: %d\n", flags[flagIndex], tdPathname, numRows);
        FILE* tdFile = fopen(tdPathname, "r"); // open file at path for reading
        printf("    Opening %s for reading...\n", tdPathname);
        *(td_ptr + numFlags) = readCSVData(numRows, max_layers, tdFile, &flippedFlag, path_flags_found); // collect
        // Next cycle
        (void)fgets(line,sizeof(line),in_file); // get next line
        numFlags++; // increment number of flags found
    }
    // TO DO: this for-loop must only collect time-independent data
    // use *path_flags_found...
    for (layer=0; layer<max_layers; layer++) {
        if (layer != 0) { (void)fgets(line,sizeof(line),in_file); } // don't read line on first scan, already done

        long double ldTemp[5]; // array of 5 temporary long doubles
        long double* matprop_ptrs[5] = { &layers[layer].rad, // array of pointers to various matprop per layer
                                        &layers[layer].dens,
                                        &layers[layer].shear,
                                        &layers[layer].bulk,
                                        &layers[layer].visc };

        // Create empty string and fill based on which variables are time-independent
        char stitchedLine[25] = "";
        stitchLine(stitchedLine, path_flags_found);

        // Read in detected data into temporary variables
        (void)sscanf(line,stitchedLine,ldTemp,ldTemp+1,ldTemp+2,ldTemp+3,ldTemp+4);

        // Conditional matprop assignments
        int j;
        int k = 0;
        for (j = 0; j < 5; j++) {
            *(matprop_ptrs[j]) = !(*path_flags_found & (1 << j)) ? ldTemp[k++] : TINY;
            // If time-independent, store detected data
            // If time-dependent, set to small value (~0); will average out later
        }
        
        /*for (layer=0; layer<max_layers; layer++) {
        (void)fgets(line,sizeof(line),in_file);
        (void)sscanf(line,"%Le %Lf %Le %Le %Le",&temp1,&temp2,&temp3,&temp4,&temp5);
        layers[layer].rad = temp1;
        layers[layer].dens = temp2;
        layers[layer].shear = temp3;
        layers[layer].bulk = temp4;
        layers[layer].visc = temp5; 
        }*/
    }
}

/************************************************************
 * set_grids -- Sets up the grids and derived quantities    *
 *		for each layer.				                        *
 *							                                *
 * Parameters						                        *
 *	max_layers	-- the total number of layers	            *
 ************************************************************/

void set_grids(int max_layers, char* rheology_flag) {

    int layer;	    /* the current layer */
    int i,j,k;	    /* counters */

    /* Allocate memory */
    
    grids.theta = (long double *)malloc((grids.num_theta+1)*sizeof(long double));
    grids.phi = (long double *)malloc((grids.num_phi+1)*sizeof(long double));
    grids.radius = (long double *)malloc((grids.num_r+1)*sizeof(long double));
    grids.etheta = (long double *)malloc((grids.num_theta+1)*sizeof(long double));
    grids.ephi = (long double *)malloc((grids.num_phi+1)*sizeof(long double));
    grids.eradius = (long double *)malloc((grids.num_r+1)*sizeof(long double));
    grids.ice_radius = (long double *)malloc((grids.num_r+1)*sizeof(long double));

    /* Initialize theta, phi, radius vectors */
    fprintf(stdout,"%d x %d x %d nodes\n",grids.num_theta,grids.num_phi,grids.num_r);
    grids.dth = M_PI / (grids.num_theta - 1);
    for(i=0;i<grids.num_theta;i++) {
	grids.theta[i] = i*grids.dth;
	if(i>=1)
	    grids.etheta[i] = 0.5 * (grids.theta[i-1] + grids.theta[i]);
    }
    /* Add endpoints to element vector */
    grids.etheta[0] = grids.theta[0];
    grids.etheta[grids.num_theta] = grids.theta[grids.num_theta-1];

    grids.df = 2.0*M_PI / (grids.num_phi - 1);
    for(j=0;j<grids.num_phi;j++) {
	grids.phi[j] = j*grids.df;
	if(j>=1)
	    grids.ephi[j] = 0.5 * (grids.phi[j-1] + grids.phi[j]);
    }
    /* Add endpoints to element vector */
    grids.ephi[0] = grids.phi[0];
    grids.ephi[grids.num_phi] = grids.phi[grids.num_phi-1];

    /* Take radius vector as radius of each layer.
     * For element grid, eradius[1] = (radius[2] + radius[3] / 2)
     */
    for(k=0;k<grids.num_r;k++) {
	grids.radius[k] = layers[k].rad;
	/* Set element r grid, start with 1 at bottom, max_layers-3 at top */
	/* Ice r grid has two fewer elements, no silicate or ocean layer */
	if (k>=1) {
	    grids.eradius[k] = 0.5 * (grids.radius[k-1] + grids.radius[k]);
//	    (void)fprintf(stderr,"element %d eradius %Le ",k,grids.eradius[k]);
	    if (k>=3) {
		grids.ice_radius[k-2] = grids.eradius[k];
//		(void)fprintf(stderr,"ice_radius %Le",grids.ice_radius[k-2]);
		}
//	    (void)fprintf(stderr,"\n");
	}
//	(void)fprintf(stderr,"node %d radius %Le\n",k,grids.radius[k]);
    }

    /* For convenience add surface to end of eradius and ice_radius arrays */
    grids.eradius[max_layers] = layers[max_layers-1].rad;
    grids.eradius[0] = layers[0].rad;
    grids.ice_radius[max_layers-2] = layers[max_layers-1].rad;
    grids.ice_radius[0] = layers[2].rad;
    grids.r_outer = layers[max_layers-1].rad;
    grids.r_inner = layers[2].rad;
    grids.dr = (grids.r_outer - grids.r_inner) / (grids.num_r - 3);

    /* Compute derived complex and gravitational quantities */
    for (layer=0; layer<max_layers; layer++) {
	/* Viscoelastic "rigidity" */
	if (layer==0) {
        layers[layer].rigidity = 0.0;
    }
    else { // if (*rheology_flag & ANDRADE_RHEOLOGY)
        long double omega_to_alpha_pow = pow(planet.freq, andrade.alpha);
        layers[layer].rigidity = ((1 / layers[layer].shear + andrade.beta * andrade.C / omega_to_alpha_pow) + \
                                  I * (1 / (layers[layer].visc*planet.freq) + andrade.beta*andrade.S / omega_to_alpha_pow)) \
                                 / ((1 / layers[layer].shear + andrade.beta * andrade.C / omega_to_alpha_pow)*(1 / layers[layer].shear + andrade.beta * andrade.C / omega_to_alpha_pow) + \
                                 (1 / (layers[layer].visc*planet.freq) + andrade.beta*andrade.S / omega_to_alpha_pow)*(1 / (layers[layer].visc*planet.freq) + andrade.beta*andrade.S / omega_to_alpha_pow));
    }
	/*else {
	    layers[layer].rigidity = 
		(layers[layer].shear * planet.freq*planet.freq * 
		    layers[layer].visc*layers[layer].visc) / 
		(layers[layer].shear*layers[layer].shear 
		    + planet.freq*planet.freq 
			* layers[layer].visc*layers[layer].visc)
	    + I *
		(layers[layer].shear*layers[layer].shear * planet.freq 
		    * layers[layer].visc) /
		(layers[layer].shear*layers[layer].shear 
		    + planet.freq*planet.freq 
			* layers[layer].visc*layers[layer].visc) ;
    }*/
	//layers[layer].rigidity = layers[layer].shear;
	/* Lame Parameter */
	layers[layer].lambda = layers[layer].bulk 
				- (2.0/3.0)*layers[layer].rigidity;
	/* Grav. Parameter */
	layers[layer].A = (4.0/3.0) * PI * G * layers[layer].dens;

	/* Calc. gravity at top of each layer */
	layers[layer].g = layers[layer].A * layers[layer].rad;
	/* Add contribution from each layer beneath */
	for(i = (layer-1); i >= 0; i--) { 
	    layers[layer].g += (layers[i].A - layers[i+1].A) 
		* powf(layers[i].rad,3) / powf(layers[layer].rad,2);
	}
	
/*	(void)fprintf(stderr,"%Lg + %Lgi\t%Lg + %Lgi\t%Lg\t%Lg\n",
		creall(layers[layer].rigidity),cimagl(layers[layer].rigidity),
		creall(layers[layer].lambda),cimagl(layers[layer].lambda),
		layers[layer].A,layers[layer].g);
*/	
		
    }

    /* Compute volume */
    planet.vol = (4.0*PI/3.0) * pow(grids.r_outer,3);

/*    (void)fprintf(stderr,"%s\n",dir);
    strcpy(directory,dir);
    directory[strlen(directory)-1] = '\0';*/
    /*(void)fprintf(stderr,"%s\n",directory);*/
}


void print_stuff(char directory[]) {

    int i,j,k;
    char filename[80];			/* Test File name */

    FILE *out_file1;

    /* check outputs */
    (void)fprintf(stdout,"%d x %d x %d nodes\n",
	    grids.num_theta,grids.num_phi,grids.num_r);
    (void)fprintf(stderr,"dtheta = %.1Lf deg, dphi = %.1Lf deg , dr = %.1Lf km\n",
	    grids.dth*DEGRAD,grids.df*DEGRAD,grids.dr/1000.0);
    (void)fprintf(stderr,"r_inner %.2Lf km, r_outer %.2Lf km\n",
	    grids.r_inner/1000.0,grids.r_outer/1000.0);
    (void)fprintf(stderr,"eta %.1Le Pa s, Ea %.1Lf kJ/mol, k %.1Lf W/(m K)\n",
	    planet.visc0,planet.E_a/1000.0,planet.k);
    (void)fprintf(stderr,"freq %.3Le s-1, e %.4Lf\n, i %.4Lf\n",
		planet.freq,planet.ecc,planet.ob);
    (void)fprintf(stderr,"Tsurf %.1Lf K, Tbot %.1Lf K, Fbot %.1Lf mW/m^2\n",
	    bc.surftemp,bc.bottemp,bc.flux*1000.0);

    /* Print coordinates */
    sprintf(filename,"%s/coords.dat",directory);
    out_file1 = fopen(filename,"w");
    for(i=0;i<=grids.num_theta-1;i++) 
	for(j=0;j<=grids.num_phi-1;j++)
	    for(k=0;k<=grids.num_r-1;k++) 
		(void)fprintf(out_file1,"%Lf %Lf %Lf\n",
			grids.theta[i],grids.phi[j],grids.radius[k]);
    fclose(out_file1);

    sprintf(filename,"%s/coords_elem.dat",directory);
    out_file1 = fopen(filename,"w");
    for(i=0;i<=grids.num_theta;i++)
	for(j=0;j<=grids.num_phi;j++)
	    for(k=0;k<=grids.num_r;k++)
		(void)fprintf(out_file1,"%Lf %Lf %Lf\n",
			grids.etheta[i],grids.ephi[j],grids.eradius[k]);
    fclose(out_file1);

    sprintf(filename,"%s/coords_ice_elem.dat",directory);
    out_file1 = fopen(filename,"w");
    for(i=0;i<=grids.num_theta;i++)
	for(j=0;j<=grids.num_phi;j++)
	    for(k=0;k<=(grids.num_r-2);k++)
		(void)fprintf(out_file1,"%Lf %Lf %Lf\n",
			grids.etheta[i],grids.ephi[j],grids.ice_radius[k]);
    fclose(out_file1);

}

/************************************************************
 * compute_visc -- get the temperature-dependent viscosity  *
 *		finite-difference CS scheme and calculate           *
 *		the misfit at the surface.		                    *
 *							                                *
 * Parameters						                        *
 *	max_layers  -- Number of layers in planet	            *
 *	T	    -- Temperature			                        *
 *	directory   -- data directory			                *
 *	tc_iter	    -- which iteration is this		            *
 ************************************************************/

void compute_visc(int max_layers, long double **T, char directory[], 
	int tc_iter) {

    FILE *out_file;	    /* Output file */
    char output_file[80];   /* Full path of heating_file */
    int i,k,layer;	    /* Counters */
    int nr,nth;	    
    long double th,r;
    long double dth,dr;	    /* deltas */
    long double thi, rk;    /* current coords */
    long double Tdep;	    /* Temperature-dependence */
    long double temp;	    /* Horiz avg. Temperature */
    long double visc;	    /* Horiz avg. Viscosity */
    long double heat;	    /* Horiz avg. heating */
    long double vol;

    nr = grids.num_r - 2;
    nth = grids.num_theta;

    dth = grids.dth;
    dr = grids.dr;

    sprintf(output_file,"%s/radials.%d",directory,tc_iter);
    out_file = fopen(output_file,"w");
    fprintf(stderr,"%s\n",output_file);


    for (k=0;k<=nr;k++) {
	r = grids.ice_radius[k];
	/* Calculate horizontally averaged temperature, heating */
	temp = 0.0;
	heat = 0.0;
	vol = 0.0;
	for (i=1;i<=(nth-1);i++) {
	    th = grids.etheta[i];
	    temp += T[i][k] * sin(th);
	    heat += heating_axi_ice[i][k] * sin(th);
	    vol += sin(th);
//	    if (k==1)
//		(void)fprintf(stderr,"T %d %Lf\n",i,T[i][k]);
	}
//	if (k==1)
//	    (void)fprintf(stderr,"T %Lf %Lf %Lf\n",temp/vol,vol,dth);
	/* Normalize */
	/*temp *= (0.5 * dth);  
	vol *= (0.5 * dth);
	heat *= (0.5 * dth);*/

	temp /= vol;
	heat /= vol;

	/* Calculate viscosity at depth k */
	Tdep = (planet.E_a / R_G) * (1.0/temp - 1.0/bc.bottemp);
	visc = planet.visc0 * exp(Tdep);
//	(void)fprintf(stderr,"%d %Lf %Lf %Lf %Le\n",k,vol,temp,Tdep,visc);
	
	/* Set viscosity of a layer to horizontal average computed here */
	/* Ignore boundaries */
	if (k>=1 && k<=(max_layers-1))
	    /* Shift up two layers, leave sil, ocean alone */
	    layers[k+2].visc = visc;

	(void)fprintf(out_file,"%Le %Le %Le %Le\n",r,temp,heat,visc);
    // TODO: store this!
    }

    fclose(out_file);
    /* Verify new viscosities */
    for (k=0;k<=(max_layers-1);k++)
//	(void)fprintf(stderr,"layer %d visc %Le Pa s\n",k,layers[k].visc);


    /* Update derived complex and gravitational quantities */
    /* Layers 2 and up (Ocean and Ice ) */
    for (layer=1; layer<=(max_layers-1); layer++) {
	layers[layer].rad = grids.radius[layer];

	/* Viscoelastic "rigidity" */
	layers[layer].rigidity = 
		(layers[layer].shear * planet.freq*planet.freq * 
		    layers[layer].visc*layers[layer].visc) / 
		(layers[layer].shear*layers[layer].shear 
		    + planet.freq*planet.freq 
			* layers[layer].visc*layers[layer].visc)
	    + I *
		(layers[layer].shear*layers[layer].shear * planet.freq 
		    * layers[layer].visc) /
		(layers[layer].shear*layers[layer].shear 
		    + planet.freq*planet.freq 
			* layers[layer].visc*layers[layer].visc) ;
	/* Lame Parameter */
	layers[layer].lambda = layers[layer].bulk 
				- (2.0/3.0)*layers[layer].rigidity;

	/* Calc. gravity at top of each layer */
	layers[layer].g = layers[layer].A * layers[layer].rad;
	/* Add contribution from each layer beneath */
	for(i = (layer-1); i >= 0; i--) { 
	    layers[layer].g += (layers[i].A - layers[i+1].A) 
		* powf(layers[i].rad,3) / powf(layers[layer].rad,2);
	}
	/*
	(void)fprintf(stderr,"%Lg + %Lgi\t%Lg + %Lgi\t%Lg\t%Lg\n",
		creall(layers[layer].rigidity),cimagl(layers[layer].rigidity),
		creall(layers[layer].lambda),cimagl(layers[layer].lambda),
		layers[layer].A,layers[layer].g);
	*/
		
    }

/*    for (layer=0; layer<=(max_layers-1); layer++)
	(void)fprintf(stderr,"layer %d, r %Le visc %Le m\n",layer,layers[layer].rad,layers[layer].visc);
*/
}
