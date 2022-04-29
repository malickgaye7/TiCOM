/* in this file define the contents of the VISC_OPT data structure
   which is used to store information used to create predefined 
   viscosity fields, those determined from prior input, those
   related to temperature/pressure/stress/anything else. */


struct VISC_OPT {
    void (* update_viscosity)();
  
    int update_allowed;		/* determines whether visc field can evolve */
    int EQUIVDD;			/* Whatever the structure, average in the end */
    int equivddopt;
    int proflocx;			/* use depth dependence from given x,y location */
    int proflocy;
    int SMOOTH;
    int smooth_cycles;
  

    char STRUCTURE[20];		/* which option to determine viscosity field, one of .... */
    int FROM_SYSTEM;
    int FROM_FILE;
    int FROM_SPECS;
  
				/* System ... */
    int RHEOL;			/* 1,2 */
    int rheol_layers;
    int num_mat;

    int nlm;
    int n410;
    int nlith;
    int ncrust1;
    int ncrust2;
    int ndd;
    float zlm;
    float z410;
    float zlith;
    float zcrust1;
    float zcrust2;
    float zdd;

    int FREEZE;
    float freeze_thresh;
    float freeze_value;

    int MAX;
    float max_value;
    int MIN;
    float min_value;

    int SDEPV;
    float sdepv_misfit;
    int sdepv_normalize;
    float sdepv_expt[40];
    float sdepv_trns[40];

    int TDEPV;
    int TDEPV_AVE;
    float N0[40];
    float E[40],T0[40];
    float T[40],Z[40];

    int weak_blobs;
    float weak_blobx[40];
    float weak_bloby[40];
    float weak_blobz[40];
    float weak_blobwidth[40];
    float weak_blobmag[40];
   
    int weak_zones;
    float weak_zonex1[40];
    float weak_zoney1[40];
    float weak_zonez1[40];
    float weak_zonex2[40];
    float weak_zoney2[40];
    float weak_zonez2[40];
  
    float weak_zonewidth[40];
    float weak_zonemag[40];
  
    int guess;
    char old_file[100];
				/* Specification info */
  
				/* Prespecified viscosity parameters */
    char VISC_OPT[20];

    int layers;			/* number of layers with properties .... */
    float layer_depth[40];
    float layer_visc[40];

    int SLABLVZ;			/* slab structure imposed on top of 3 layer structure */
    int slvzd1,slvzd2,slvzd3;	        /* layer thicknesses (nodes) */
    int slvzD1,slvzD2;		        /* slab posn & length */
    float slvzn1,slvzn2,slvzn3,slvzN;   /* viscosities */

    int COSX;
    float cosx_epsilon;
    float cosx_k;
    int cosx_exp;
 
    int EXPX;
    float expx_epsilon;
 
    /* MODULE BASED VISCOSITY VARIATIONS */

    int RESDEPV;
    float RESeta0[40];

    int CHEMDEPV;
    float CH0[40];
    float CHEMeta0[40];
  
} viscosity;
