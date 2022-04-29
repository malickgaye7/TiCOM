 	/* This file contains the definitions of variables which are passed as arguments */
	/* to functions across the whole filespace of CITCOM. #include this file everywhere !*/

#include <assert.h>
#include <stdio.h>
#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif
#ifndef COMPLEX_H
#define COMPLEX_H
#include <complex.h>
#endif

#if defined(__osf__) 
void *Malloc1();
#endif

#define Malloc0(a) Malloc1((a),__FILE__,__LINE__)


/* #define Malloc0 malloc */
#define LIDN 0x1
#define VBX 0x2
#define VBZ 0x4
#define VBY 0x8
#define TBX 0x10
#define TBZ 0x20
#define TBY 0x40
#define TZEDGE 0x80
#define TXEDGE 0x100
#define TYEDGE 0x200
#define VXEDGE 0x400
#define VZEDGE 0x800
#define VYEDGE 0x1000
#define INTX 0x2000
#define INTZ 0x4000
#define INTY 0x8000
#define SBX 0x10000
#define SBZ 0x20000
#define SBY 0x40000
#define FBX 0x80000
#define FBZ 0x100000
#define FBY 0x200000

#define CBX 0x400000
#define CBZ 0x800000
#define CBY 0x1000000
#define HBX 0x2000000
#define HBZ 0x4000000
#define HBY 0x8000000

#define OFFSIDE 0x10000000

#define SKIP 0x1000000
#define SKIPID 0x1
#define ZEROID 0x2

#define REFINE1 0x1
#define REFINE2 0x2
#define GREFINE1 0x4
#define GREFINE2 0x8

#define LIDE 1

#ifndef COMPRESS_BINARY
#define COMPRESS_BINARY "/usr/bin/compress"
#endif

#define MAX_LEVELS 12
#define MAX_F    10
#define MAX_S    30

/* Macros */

#define max(A,B) (((A) > (B)) ? (A) : (B))
#define min(A,B) (((A) < (B)) ? (A) : (B))
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define QPI 3.1415926535897932384626433832795

typedef float higher_precision;  /* matrix coeffs etc */
typedef double higher_precision1; /* intermediate calculations for finding above coeffs */


/* Common structures */

struct Segment  {

                int zlayers;
                int nzlayer[40];
                float zzlayer[40];
                int xlayers;
                int nxlayer[40];
                float xxlayer[40];
                int ylayers;
                int nylayer[40];
                float yylayer[40];

};


struct Rect {
    int numb;
    char overlay[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float y1[40];
    float y2[40];
    float halfw[40];
    float mag[40];
} ;
 

struct Circ {
    int numb;
    char overlay[40];
    float x[40];
    float z[40];
    float y[40];
    float rad[40];
    float mag[40];
    float halfw[40];
};


struct Harm {
    int numb;
    int harms;
    char overlay[40];
    float off[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float y1[40];
    float y2[40];
    float kx[20][40];
    float kz[20][40];
    float ky[20][40];
    float ka[20][40];
    float phx[20][40];
    float phz[20][40];
    float phy[20][40];
};

struct Erfc {
 int numb;


};

struct RectBc {
    int numb;
    char norm[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float halfw[40];
    float mag[40];
} ;
 

struct CircBc {
    int numb;
    char norm[40];
    float intercept[40];
    float x[40];
    float z[40];
    float rad[40];
    float mag[40];
    float halfw[40];
};


struct PolyBc {
    int numb;
    int order;
    char norm[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float ax[20][40];
    float az[20][40];
};


struct HarmBc {
    int numb;
    int harms;
    char norm[40];
    float off[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float kx[20][40];
    float kz[20][40];
    float ka[20][40];
    float phx[20][40];
    float phz[20][40];
 };


struct Shape_function_dA  {
  double vpt[8];
  double spt[4];
  double ppt[1]; };

struct Shape_function1_dA  {
  double vpt[6*4];
  double ppt[6*1]; };

struct Shape_function1 	{ 
    double vpt[4*4];  /* node & gauss pt */
    double ppt[4*1];  };

struct Shape_function 	{ 
    double vpt[8*8];  /* node & gauss pt */
    double spt[8*4];  /* node & gauss pt */
    double ppt[8*1];  };

struct Shape_function_dx 	{ 
    double vpt[3*8*8]; /* dirn & node & gauss pt */
    double spt[3*8*4]; /* dirn & node & gauss pt */
    double ppt[3*8*1];  };

struct Shape_function1_dx 	{ 
    double vpt[2*4*4]; /* dirn & node & gauss pt */
    double ppt[2*4*1];  };

struct EG { 
    higher_precision g[24][1]; };

struct EK2 { 
    double k[8*8]; };

struct EK { 
    double k[24*24]; };

struct MEK { 
    double nint[9]; };
 
struct NK {
    higher_precision *k;
    int *map;
};

struct COORD {
    float centre[4];
    float size[4];
    float recip_size[4];
    float area;   } ;

struct SUBEL { 
    int sub[9];   };
			
struct ID  { 
    int doff[6];	}; /* can  be 1 or 2 or 3 */
struct IEN {
    int node[9];	};
struct FNODE {
    float node[9];	};
struct SIEN {
    int node[5];	};
struct LM  { 
    struct { int doff[4]; } node[9]; } ;

struct NEI {
    int *nels;
    int *lnode;
    int *element; };

struct BOUND  { 
    int bound[8];	}; 

struct IBM_DX   {
    float *x1;
    float *x2;
    int nox;
    int noz;
    };

struct MESH_DATA {/* general information concerning the fe mesh */ 
    int nsd;        /* Spatial extent 1,2,3d*/
    int dof;        /* degrees of freedom per node */
    int levmax;     
    int levmin;
    int levels;
    int mgunitx;
    int mgunitz;
    int mgunity;
    int NEQ[MAX_LEVELS];	/* All other values refer to the biggest mesh (& lid)  */
    int NNO[MAX_LEVELS];
    int NNOV[MAX_LEVELS];
    int NLNO[MAX_LEVELS];
    int NPNO[MAX_LEVELS];
    int NEL[MAX_LEVELS];
    int NOX[MAX_LEVELS];
    int NOZ[MAX_LEVELS];
    int NOY[MAX_LEVELS];
    int NNX[MAX_LEVELS][4];
    int ELX[MAX_LEVELS];
    int ELZ[MAX_LEVELS];
    int ELY[MAX_LEVELS];
    int LNDS[MAX_LEVELS];
    int LELS[MAX_LEVELS];
    int neqd;
    int neq;
    int nno;
    int nnov;
    int nlno;
    int npno;
    int nel;
    int snel;
    int elz;
    int ely;
    int nnx[4]; /* general form of ... */
    int nox;
    int elx;
    int noz;
    int noy;
    int *exs;
    int ezs;
    int eys;
    int nxs;
    int nzs;
    int nys;
    int nmx;
    int nsf; /* nodes for surface observables */
    int esf; /* elements for surface observables */
    int toptbc,topcbc;
    int bottbc,botcbc;
    int topvbc;
    int botvbc;
    int sidevbc;

    char topvbc_file[100];
    char botvbc_file[100];
    char sidevbc_file[100];
    char gridfile[4][100];


    int periodic_x;
    int periodic_y;
    float layer[4];			/* dimensionless dimensions */
    float lidz;
    float bl1width[4],bl2width[4],bl1mag[4],bl2mag[4];
    float hwidth[4],magnitude[4],offset[4],width[4]; /* grid compression information */ 
    int fnodal_malloc_size;
    int dnodal_malloc_size;
    int feqn_malloc_size;
    int deqn_malloc_size;
    int bandwidth;
    int null_source;
    int null_sink;
    int matrix_size[MAX_LEVELS];

} ;

struct HAVE {    /* horizontal averages */
    float *T;
    float *Vi;
    float *r_Ri;
    float *i_Ri;
    float *Tadi;
    float *Rho;
    float *f;
    float *F;
    float *vrms;
    float *V[4];
    float *Tprev;
};

struct SLICE {    /* horizontally sliced data, including topography */
    float *tpg;
    float *tpgb;
    float *grv;
    float *geo;
    float *geok;
    float *grvk;
    float *grvb;
    float *geob;
    float *geobk;
    float *grvbk;
    float *tpgk;
    float *tpgbk;
    float *shflux;
    float *bhflux;
    float *cen_mflux;
    float *vxsurf[3];    /* surface velocity vectors */
    float *vline;        /* for kernels, velocity at force term */
    float *vlinek;
    float *tpglong;
    float *tpgblong;
    float *grvlong;
    float *geolong;
    float *grvblong;
    float *geoblong;
		float *melt;
		float *new_melt;
		float *impact_melt;

    int minlong;
    int maxlong;
  };

struct IMPACTS {
    char heating_file[100]; /* file containing impact heating */

    int number;
    int flag[40];       /* flag for whether heating has been applied */
    int now;       /* flag for whether impact is occuring now */
    int heat_from_file; /* flag to read in impact heating from file */

    float t[40];        /* time of impacts */
    float th[40];       /* co-lat and long. */
    float f[40];
    float size[40];     /* size of heated region */
    float dT;       /* Temp increase at center of heated region */
    float v;            /* Impact velocity */
    float *H_t;       /* transient depth */
};



struct SPHERE  {
    float ro_dim;
    float ro;
    float ri;
    float rcomp;
    float rcore;
    float deltarb;  /* Change in shell thickness */
    };

struct BAVE {
    float T;
    float Vi;
    double V[4]; };


struct TOTAL {
    float melt_prod; /* Total amount of melt produced */
    float vol;       /* mantle volume */
    float bulk_comp; /* average composition of mantle */
    float bulk_comp_prev; /* average composition of mantle */  };

struct MONITOR {
    char node_output[100][6];  /* recording the format of the output data */
    char sobs_output[100][6];  /* recording the format of the output data */
    int node_output_cols;
    int sobs_output_cols;

    int solution_cycles;

    float  time_scale;
    float  length_scale;
    float  viscosity_scale;
    float  geoscale;
    float  tpgscale;
    float  grvscale;
  
    float  delta_v_last_soln;
    float  elapsed_time;
    float  elapsed_time_vsoln;
    float  elapsed_time_vsoln1;
    float  reference_stress;
    float  incompressibility;
    float  vdotv;
    float  nond_av_heat_fl;  
    float  nond_av_adv_hfl;  
    float  cpu_time_elapsed;
    float  cpu_time_on_vp_it;
    float  cpu_time_on_forces;
    float  cpu_time_on_mg_maps;
    float  tpgkmag;
    float  grvkmag;
   
    float  Nusselt;
    float  Vmax;
    float  Vsrms;
    float  Vrms;
    float  Vrms_surface;
    float  Vrms_base;
    float  F_surface;
    float  F_base;
    float  Frat_surface;
    float  Frat_base;
    float  T_interior;
    float  T_maxvaried;
    float  Sigma_max;
    float  Sigma_interior;
    float  Vi_average;

    float deltah;    /* total change in shell thickness */ 
};

struct CONTROL {
    int PID;

    char output_written_external_command[500];   /* a unix command to run when output files have been created */

    int ORTHO,ORTHOZ;   /* indicates levels of mesh symmetry */
    char B_is_good[MAX_LEVELS];  /* general information controlling program flow */
    char Ahat_is_good[MAX_LEVELS];  /* general information controlling program flow */
    char old_P_file[100];
    char data_file[100];
    char data_file1[100];

    char which_data_files[1000];
    char which_horiz_averages[1000];
    char which_running_data[1000];
    char which_observable_data[1000];
  
    char PROBLEM_TYPE[20]; /* one of ... */
    int KERNEL;
    int stokes;
    int CONVECTION;
    char comp_adv_method[20];
    int SLAB;	
    char GEOMETRY[20]; /* one of ... */
    int CART2D;
    int CART2pt5D;
    int CART3D;
    int AXI;	 
    char SOLVER_TYPE[20]; /* one of ... */
    int DIRECT;
    int CONJ_GRAD;
    int NMULTIGRID;
    int EMULTIGRID;
    int DIRECTII;
    char NODE_SPACING[20]; /* turns into ... */
    int GRID_TYPE;
    int COMPRESS;
    int DX;
    int CONMAN;
    int visc_heating;
    int adi_heating;
    int latent_heating;
    int tidal_heating;
    int shear_heating;
    int despin;
		int despun;
 
    int composition;
		int melting;
    int freezing;
    float comp_diff;
    float z_comp;
    float Q0ER;
 
    int surf_temp_var; 
    int secular; 

    int dfact;
    double penalty;
    int augmented_Lagr;
    double augmented;
    int macroele;
    int faults;
    int NASSEMBLE;
    int comparison;
    int crust;
    int restart;
    int restart_frame;
    float plate_vel;

    float tole_comp;
  
    float sob_tolerance;
 
    double Ra_temp,Ra_comp,Ra_comp_a; 
    float Ra_670,clapeyron670,transT670,width670; 
    float Ra_410,clapeyron410,transT410,width410; 
    float Ts; 
    float VBXtopval;
    float VBXbotval;
    float VBYtopval;
    float VBYbotval;

    float TBCtopval,CBCtopval;
    float TBCbotval,CBCbotval;

    float Q0;
    float Qc;
    float despin_timescale_pred;
    float despin_timescale;

		float restart_age;
    
    int dimensionalize;
    int precondition;
    int vprecondition;
    int keep_going;
    int v_steps_low;
    int v_steps_high;
    int v_steps_upper;
    int max_vel_iterations;
    int p_iterations;
    int max_same_visc;
    float max_res_red_each_p_mg;
    float sub_stepping_factor;
    int mg_cycle;
    int true_vcycle;
    int down_heavy;
    int up_heavy;
    int depth_dominated;
    int eqn_viscosity;
    int eqn_zigzag;
    int verbose;
    double accuracy;
    double vaccuracy; 
   
    int total_iteration_cycles;
    int total_v_solver_calls;
    
    int record_every;
    int record_all_until;

    int print_convergence;
    int sdepv_print_convergence;

     /* modules */
    int MELTING_MODULE;
    int CHEMISTRY_MODULE;
};

struct Radioactive_Heat {
    int num;
    double total;
    double concen_u;
    double percent[10];
    double heat_g[10];
    double decay_t[10];
    double concen[10];
};

struct DATA {  
    float  layer_km;
    float   grav_acc;
    float   therm_exp;
    float   therm_exp_factor;
    float   visc_factor;
    float   Cp;
    float   Hf;

    float  starting_time;
    float  crust_rad_enhanced;
    float  mantle_rad;
    float  crust_rad;
    float  crust_mantle_vol_ratio;

    float  disptn_number;
    float  therm_diff;
    float  therm_diff_factor;
    float  therm_cond;
    float   density;
    float   density_lm;
    float  res_density;
    float  res_density_X;
    float   melt_density;
    float   density_above;
    float   density_core;
    float   gas_const;
    float   surf_heat_flux;
    double  ref_viscosity;
    float   melt_viscosity;
    float   permeability;
    float   grav_const;
    float  surf_temp;
    float   youngs_mod; 
    float   Te;
    float   ref_temperature;
    float   DeltaT;
    float   Tsurf;
    float   T_sol0;
    float   T_adi0;
    float   T_adi1;

    double  rigidity;
    double  frequency;
    double  rot_init;
    double  rot_final;
		double  rotation;
		double	flattening;
    double  mass_primary;
    double  mass;
    long double  semimajor_axis;
		long double	mean_motion;
		double	moi;
		double	despin_rate;
		double	Ediss;
		double	Ediss_total;
		
    float   delta_S;
    float   dTsol_dz;
    float   dTsol_dF;
    float   dT_dz; };

struct NEWVARS {
    double *P,*U;
    float *XX[MAX_LEVELS];
    float *oldX;
    float *oldecocentre;
    float *XP;
    float *T,*C;
    float *Fm;
    float *heating_tidal, *tidal_visc;
    float *solidus;
    float *Vi,*EVi;
    int *mat;	        /* properties of mat */
    int coords;       /* flag to output coords again */
};

struct All_variables {     
#ifndef __CONVECTION_VARIABLES_H__
#define __CONVECTION_VARIABLES_H__
#include "Convection_variables.h"
#endif

#ifndef __VISCOSITY_DESCRIPTIONS_H__
#define __VISCOSITY_DESCRIPTIONS_H__
#include "viscosity_descriptions.h"
#endif

#ifndef __TEMPERATURE_DESCRIPTIONS_H__
#define __TEMPERATURE_DESCRIPTIONS_H__
#include "temperature_descriptions.h"
#endif

#ifndef __ADVECTION_H__
#define __ADVECTION_H__
#include "advection.h"
#endif

    FILE *fp;
    FILE *fpdebug;
    FILE *fpq;
    FILE *filed[20];
    struct HAVE Have;
    struct BAVE Bulkave;
    struct TOTAL Total;
    struct MESH_DATA mesh;
    struct CONTROL control;
    struct MONITOR monitor;
    struct DATA data;
    struct SLICE slice;
    struct COORD *eco;
    struct IBM_DX ibm_dx;
    struct IEN *ien;  /* global */
    struct SIEN *fault_ien;
    struct SIEN *sien;
    struct ID *id;
    struct COORD *ECO[MAX_LEVELS];
    struct IEN *IEN[MAX_LEVELS]; /* global at each level */
    struct FNODE *TWW[MAX_LEVELS]; /* */
    struct ID *ID[MAX_LEVELS];
    struct NEI NEI[MAX_LEVELS];
    struct SUBEL *EL[MAX_LEVELS];
    struct EG *elt_del[MAX_LEVELS];
    struct EK *elt_k[MAX_LEVELS];
    struct SPHERE sphere;
    struct Radioactive_Heat rad_heat;
    struct Segment segment;
		struct IMPACTS impacts;
    struct NEWVARS newvars;

    higher_precision *Eqn_k[MAX_LEVELS];  
    int *Node_map[MAX_LEVELS];
    int *Node_eqn[MAX_LEVELS];
    int *Node_k_id[MAX_LEVELS];


    double *BI[MAX_LEVELS];      /* inv of  diagonal elements of K matrix */
    double *BPI[MAX_LEVELS];
    float *V[4];                 /* velocity X[dirn][node] can save memory */
    float *V1[4];                 /* velocity X[dirn][node] can save memory */
    float *Vest[4];                 /* velocity X[dirn][node] can save memory */

    float *VO[4],*XMC[4],*XMCpred[4];
    int *C12,*CElement;
    float *C12f;

    double *P,*F,*H,*S,*U;
    float *Psi;
    float *diffusivity,*expansivity;
    float *NP;
    float *edot;               /* strain rate invariant */
    float *MASS[MAX_LEVELS],*Mass;               /* lumped mass matrix (diagonal elements) for p-g solver etc. */  
    float *tw;
    float *stress;
    float *XP[4],*X[4],*XX[MAX_LEVELS][4],*Interp[MAX_LEVELS][4];
    float *ZZ,*heatflux;
    float *T,*TE,*C,*C_prev,*CE,*CE_prev,*CE_temp,*IHeat,*buoyancy;
    float *Tdot,*Cdot;		
		float *Fm,*FmE;		
    float *heating_visc,*heating_adi,*heating_latent;		
    float *heating_tidal, *tidal_visc;
    float *heating_shear;
    double *heating_despin, power_despin_global, power_despin_ave;
    float *Fas670,*Fas410;		
    float *Fas670_b,*Fas410_b;		
		float *solidus;
		float *liquidus;
		float *lherzliq;

    float *Vi,*EVi;
    float *VI[MAX_LEVELS];	/* viscosity has to soak down to all levels */
    float *EVI[MAX_LEVELS];	/* element viscosity has to soak down to all levels */
    complex double *c_Ri,*c_ERi; /* complex rigidity */
    float *VB[4],*TB[4],*CB[4];/* boundary conditions for V,T defined everywhere */
    float *TW[MAX_LEVELS];	/* nodal weightings */

    int num_zero_resid[MAX_LEVELS];	       
    int *zero_resid[MAX_LEVELS];	       
    int *surf_element;	       
    int *surf_node;	       
    int *mat;	        /* properties of mat */
    unsigned int *node;	        /* properties of node */
    unsigned int *NODE[MAX_LEVELS];
    unsigned int *ELEMENT[MAX_LEVELS];
    unsigned int *eqn;
    unsigned int *EQN[MAX_LEVELS];
 		
    double **global_K;  /* direct solver stuff */
    double **factor_K;
    double *global_F;
    struct LM *lmd;
    struct LM *lm;
    struct LM *LMD[MAX_LEVELS];
  
    struct Shape_function1 M; /* master-element shape funtions */
    struct Shape_function1_dx Mx; 
    struct Shape_function N;
    struct Shape_function_dx Nx;
    struct Shape_function1 L; /* master-element shape funtions */
    struct Shape_function1_dx Lx; 
 
    void (* build_forcing_term)();
    void (* iterative_solver)();
    void (* next_buoyancy_field)();
    void (* obtain_gravity)();
    void (* problem_settings)();
    void (* problem_derived_values)();
    void (* problem_allocate_vars)();
    void (* problem_boundary_conds)();
    void (* problem_node_positions)();
    void (* problem_update_node_positions)();
    void (* problem_initial_fields)();
    void (* problem_update_bcs)();
    void (* special_process_new_velocity)();
    void (* special_process_new_buoyancy)();
    void (* solve_stokes_problem)(); 
    void (* solver_allocate_vars)(); 
    void (* transform)();

    float (* node_space_function[3])();
 
};
