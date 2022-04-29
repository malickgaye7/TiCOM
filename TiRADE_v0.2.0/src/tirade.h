/************************************************************
 * definitions for the Tidal Heating / Conduction	          *
 *  (tirade) package				                                *
 *							                                            *
 * There are a lot of variables that are common to many	    *
 * functions, such as the grid information and global	      *
 * physical properties.  It makes no sense to pass them	    *
 * around all day, so let's make them global to the program *
 * and define the structure types here.			                *
 *----------------------------------------------------------*
 * struct grid_info					                                *
 *	used to hold the information on the finite-	            *
 *	difference mesh spacing and co-ordinates.	              *
 *----------------------------------------------------------*
 * struct globals					                                  *
 *	used to hold information about global physical and      *
 *	orbital properties that do not change spatially.        *
 *----------------------------------------------------------*
 * struct boundaries					                              *
 *	used to hold information about boundary conditions      *
 *----------------------------------------------------------*
 * struct matprop					                                  *
 *	used to hold material proerties for a layer within      *
 *	the planet.					                                    *
 ************************************************************/

#define PI 3.1415926535	    /* Pi */
#define DEGRAD (180 / PI)   /* Radian to degree conversion */
#define TINY 1.0e-20	    /* A small number */
#define G 6.67e-11	    /* Gravitational Constant */
#define R_G 8.314	    /* Gas Constant */
#define OFFBYTE 0x00    /* Byte w/ every bit set to 0 */
#define ONBYTE ~OFFBYTE /* Byte w/ every bit set to 1 */


#define ANDRADE_ALPHA_TOO_LARGE 32.0 /* Alpha that is consered 'very large', mu goes to Maxwell solution */
// rheology_params bits
#define ANDRADE_RHEOLOGY 1<<0

typedef struct _heap_data {
    char* type; // Necessary for all types

    // ..................... Allocated when below type is specified...

    // Commonly shared variables
    long double* etheta; // heating_axi, stress_strain
    long double* eradius; // heating_axi, heating_radial, heating_radial_ice, potential, y_i, sensparam, stress_strain
    long double* heating_axi; // heating_axi
    long double* heating_radial; // heating_radial, heating_radial_ice
    long double* ephi; // potential, stress_strain
    long double* potential; // potential
    long double* y_i_real; // y_i
    long double* y_i_imag; // y_i

    // Variables exclusive to love number & imk2(?)
    long double* Loveh2_real;
    long double* Loveh2_imag;
    long double* Lovek2_real;
    long double* Lovek2_imag;
    long double* phaseh2;
    long double* phasek2;

    // Variables exclusive to Tobie's sensitivity parameters
    long double* Hmu;
    long double* hHmu;
    long double* mu_real;
    long double* mu_imag;

    // Variables exclusive to stress & strain
    long double* stress_real;
    long double* stress_imag;
    long double* strain_real;
    long double* strain_imag;
    int* nodes_ordered;

} Tiradata;
// NOTE: 2D array will be collapsed to 1D when read into heap memory unless we want matrix dimensions stored...

struct unitFlags {
    // 1 is set to the corresponding bit
    // Character byte representations initialized to 00000000
    // Filled from the RIGHT
    unsigned char tUnits; // s - min - hr - day - wk - mnth - yr (7 bits)
    unsigned char radUnits; // in - ft - yrd - mi - cm - m - km (7 bits)
    unsigned char densUnits; // lb/ft^3 - g/cm^3 - kg/m^3 (3 bits)
    unsigned char shearbulkUnits; // psi - ksi - Pa - GPa (4 bits)
    unsigned char viscUnits; // lb s/ft^2 - Ns/m^2 same as Pa s and kg/ms (2 bits)
};

struct grid_info {
    int num_theta;	/* number of latitudinal points */
    int num_phi;	/* number of longitudinal points */
    int num_r;		/* number of radial points */
    long double dth;		/* Latitudinal size of element */
    long double df;		/* Longitudinal size of element */
    long double dr;		/* Radial size of element */
    long double r_inner;	/* Inner Radius */
    long double r_outer;	/* Outer Radius */
    long double *theta;	    /* Theta nodes */
    long double *phi;	    /* Phi nodes */
    long double *radius;    /* R nodes */
    long double *etheta;    /* Theta elements */
    long double *ephi;	    /* Phi elements */
    long double *eradius;   /* R elements */
    long double *ice_radius;	/* R elements for ice shell only */
};

struct globals {
    long double k;		/* Thermal Conductivity */
    long double E_a;		/* Activation Energy */
    long double visc0;	/* Reference Viscosity */
    long double vol;    /* Volume of planet */
    long double freq;	    /* Orbital frequency */
    long double ecc;	    /* Eccentricity of orbit */
    long double ob;	    /* Obliquity wrt Sun */
    long double true_fluid_cutoff_visc;
    long double true_fluid_cutoff_rigidity;
};

struct boundaries {
    long double bottemp;	/* Temperature at Bottom */
    long double surftemp;	/* Temperature at surface */
    long double flux;		/* Flux across Bottom */
};

struct matprop {	    /* Material properties for each layer */
    long double rad;	    /* Outer radius of the layer */
    long double dens;	    /* Density */
    long double shear;	    /* Elastic shear modulus */
    long double bulk;	    /* Bulk modulus */
    long double visc;	    /* Viscosity */
    long double A;	    /* Grav. parameter (4pi/3 G rho) */
    long double g;	    /* gravity at top of layer */
    long double complex rigidity;   /* Visco-elastic rigidity (computed) */
    long double complex lambda;	    /* Lame Parameter (computed) */
};

struct andrade_params {
    long double alpha;
    long double beta;
    long double C;
    long double S;
};
 
