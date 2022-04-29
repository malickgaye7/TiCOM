/* Assumes parameter list is opened and reads the things it needs. 
   Variables are initialized etc, default values are set */


#include <math.h>
//#include <malloc.h>
#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif
#include <sys/types.h>
#ifndef __ELEMENT_DEFINITIONS_H__
#define __ELEMENT_DEFINITIONS_H__
#include "element_definitions.h"
#endif

#ifndef __GLOBAL_DEFS_H__
#define __GLOBAL_DEFS_H__
#include "global_defs.h"
#endif
//#include <stdlib.h> /* for "system" command */
#include <strings.h>
#include <complex.h>

void set_convection_defaults(E)
     struct All_variables *E;
{
    void PG_timestep_with_melting();
    void PG_timestep();
    void PG_timestep_particle();
    void read_convection_settings();
    void convection_derived_values();
    void convection_allocate_memory();
    void convection_boundary_conditions();
    void node_locations();
    void convection_initial_fields();
    void twiddle_thumbs();

    int input_int(char*, int*, char*);
    int input_string(char*, char*, char*);
    int input_float(char*, float*, char*);

   input_int("composition",&(E->control.composition),"0");
   input_int("melting",&(E->control.melting),"0");
   input_float("comp_diffusivity",&(E->control.comp_diff),"0");
   input_string("comp_adv_method",E->control.comp_adv_method,NULL);


    if (E->control.composition && strcmp(E->control.comp_adv_method,"field")==0)
      E->next_buoyancy_field = PG_timestep;
    else if (E->control.composition && strcmp(E->control.comp_adv_method,"particle")==0)  {
      E->next_buoyancy_field = PG_timestep_particle;
      }
    else {
      E->next_buoyancy_field = PG_timestep;
      }

 
    E->special_process_new_buoyancy = twiddle_thumbs; 
    E->problem_settings = read_convection_settings;
    E->problem_derived_values = convection_derived_values;
    E->problem_allocate_vars = convection_allocate_memory;
    E->problem_boundary_conds = convection_boundary_conditions;
    E->problem_initial_fields = convection_initial_fields;
    E->problem_node_positions = node_locations;
    E->problem_update_node_positions = twiddle_thumbs;
    E->problem_update_bcs = twiddle_thumbs;

    sprintf(E->control.which_data_files,"Temp,Strf,Pres");
    sprintf(E->control.which_horiz_averages,"Temp,Visc,Vrms");
    sprintf(E->control.which_running_data,"Step,Time,");
    sprintf(E->control.which_observable_data,"Shfl");
 
return;
}

void read_convection_settings(E)
     struct All_variables *E;
    
{ 
    void advection_diffusion_parameters();
   float density_diff;
    
/* parameters */
    int input_double(char*, double*, char*);
    input_double("rayleigh",&(E->control.Ra_temp),"essential");

    E->data.ref_viscosity = E->data.grav_acc*E->data.density*E->data.therm_exp
                  *E->data.ref_temperature*E->sphere.ro_dim*E->sphere.ro_dim*E->sphere.ro_dim
                  /(E->control.Ra_temp*E->data.therm_diff);

    E->data.ref_viscosity /= (1.0 + E->data.surf_temp);

    input_double("rayleigh_comp",&(E->control.Ra_comp),"essential");

    density_diff=E->control.Ra_comp*E->data.ref_viscosity*E->data.therm_diff/(E->data.grav_acc*E->sphere.ro_dim*E->sphere.ro_dim*E->sphere.ro_dim);

    fprintf(E->fp,"Ra_temp=%.5e Ra_comp=%.5e ref_visc=%.5e density_diff=%.5e\n",E->control.Ra_temp,E->control.Ra_comp,E->data.ref_viscosity,density_diff);

    int input_boolean(char*, int*, char*);
    int input_float(char*, float*, char*);
    input_boolean("halfspace",&(E->convection.half_space_cooling),"off");
    input_float("halfspage",&(E->convection.half_space_age),"nodefault");
    
    int input_int(char*, int*, char*);
    input_int("temperature_blobs",&(E->convection.temp_blobs),"0");

    int input_float_vector(char*, int, float*);
    int input_int_vector(char*, int, int*);
    input_float_vector("temperature_blobx",E->convection.temp_blobs,E->convection.temp_blob_x);
    input_float_vector("temperature_bloby",E->convection.temp_blobs,E->convection.temp_blob_y);
    input_float_vector("temperature_blobz",E->convection.temp_blobs,E->convection.temp_blob_z);
    input_float_vector("temperature_blobsize",E->convection.temp_blobs,E->convection.temp_blob_radius);
    input_float_vector("temperature_blobDT",E->convection.temp_blobs,E->convection.temp_blob_T);
    input_float_vector("temperature_blobbg",E->convection.temp_blobs,E->convection.temp_blob_bg);
    input_int_vector("temperature_blobsticky",E->convection.temp_blobs,E->convection.temp_blob_sticky);
    
    input_int("temperature_zones",&(E->convection.temp_zones),"0");
    input_float_vector("temperature_zonex1",E->convection.temp_zones,E->convection.temp_zonex1);
    input_float_vector("temperature_zonex2",E->convection.temp_zones,E->convection.temp_zonex2);
    input_float_vector("temperature_zonez1",E->convection.temp_zones,E->convection.temp_zonez1);
    input_float_vector("temperature_zonez2",E->convection.temp_zones,E->convection.temp_zonez2);
    input_float_vector("temperature_zoney1",E->convection.temp_zones,E->convection.temp_zoney1);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zonehw",E->convection.temp_zones,E->convection.temp_zonehw);
    input_float_vector("temperature_zonemag",E->convection.temp_zones,E->convection.temp_zonemag);
    input_int_vector("temperature_zonesticky",E->convection.temp_zones,E->convection.temp_zone_sticky);
    
    input_int("num_perturbations",&(E->convection.number_of_perturbations),"0,0,32");
    input_float_vector("perturbmag",E->convection.number_of_perturbations,E->convection.perturb_mag);
    input_float_vector("ll",E->convection.number_of_perturbations,E->convection.perturb_ll);
    input_float_vector("mm",E->convection.number_of_perturbations,E->convection.perturb_mm);
    
    int input_string(char*, char*, char*);
    input_string("prevT",E->convection.old_T_file,"initialize");
    
    advection_diffusion_parameters(E);
		//fprintf(stderr,"add %d %d\n",E->advection.markers_per_ele,E->mesh.dof);
    //fprintf(stdout,"add %d %d\n",E->advection.markers_per_ele,E->mesh.dof); // DEBUG
    
    if (E->control.restart)    {
      input_int("restart_timesteps",&(E->monitor.solution_cycles),"0");
      input_string("oldfile",E->convection.old_T_file,"initialize");

			E->control.restart_frame = E->monitor.solution_cycles;
      }
		else 
			E->control.restart_frame = 0;

    return;
}

/* =================================================================
   Any setup which relates only to the convection stuff goes in here
   ================================================================= */

void convection_derived_values(E)  
     struct All_variables *E;
 
{ 

return;
}

void convection_allocate_memory(E)
     struct All_variables *E;

{ void advection_diffusion_allocate_memory();

  advection_diffusion_allocate_memory(E);

return;
}

/* ============================================ */
    
void convection_initial_fields(E)
     struct All_variables *E;

{ 
    void convection_initial_temperature();
    void convection_initial_markers();
		void define_solidus();
		void read_impacts();
    void read_tidal_heating();
    void despin();

    int report(struct All_variables *, char*);
    report(E,"convection, initial temperature");
    convection_initial_temperature(E);
		define_solidus(E);

    if(E->control.tidal_heating) {
			fprintf(stdout,"tidal\n");
			read_tidal_heating(E);
    }

    if(E->control.despin) {
			fprintf(stdout,"despin\n");
			despin(E);
    }

		if(E->impacts.number > 0) {
			fprintf(stdout,"impacts\n");
			read_impacts(E,E->impacts.number);
		}

  return; }

/* =========================================== */

void convection_boundary_conditions(E)
     struct All_variables *E;

{
    void velocity_boundary_conditions();
    void temperature_boundary_conditions();
    void temperatures_conform_bcs();
    void composition_boundary_conditions();
 
    velocity_boundary_conditions(E);      /* universal */
    temperature_boundary_conditions(E);

    temperatures_conform_bcs(E);

    composition_boundary_conditions(E);

    return;
}

/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_temperature(E)
     struct All_variables *E;
{
    int i,j,k,p,node,ii,jj;
    double temp,base,radius,radius2;
    double modified_plgndr_a(),drand48();
    FILE *fp;
    void remove_horiz_ave();
    void temperatures_conform_bcs();
    void thermal_buoyancy();
    void  process_restart_mk();
    void  process_restart_tc();
    void  convection_initial_markers();
    void  p_to_centres();
    void  reallocate_markers();
    void  get_markers_from_C();

    float a,b,c;
    float radial_temp[255];
    
    int in1,in2,in3,instance,nox,noy,noz,nfz,ok,noz2,ll,mm;
    char output_file[255], input_s[200];
    double tbase,tbase1,t1,r1,weight,para1,plate_velocity,delta_temp,age;
    double x00,x01,x02,slope,con;

    int read_previous_field();
    
    double tbl, bbl;    /* Radial Positions of Boundary Layers */
    double inttemp, adtemp;     /* Isothermal interior temperature */
    double adgrad;      /* Adiabatic Gradient */
    double T0;		/* Surf Temperature at this latitude */
	   
    const int dims=E->mesh.nsd;
    const float e_5=1.0e-5;
    
/*    sprintf(output_file,"%s/radialtemp",E->control.data_file);
    fprintf(stderr,"%s\n",output_file);
    fp=fopen(output_file,"r");
    fprintf(stderr,"opened\n");

    for(i=1;i<=E->mesh.noz;i++) {
      fprintf(stderr,"i=%d ",i);
      fgets(input_s,200,fp);
      fprintf(stderr,"i2=%d ",i);
      sscanf(input_s,"%f %f",&a,&b);
      fprintf(stderr,"i3=%d ",i);
      c = b;
      fprintf(stderr,"i4=%d ",i);
      radial_temp[i] = a;
      fprintf(stderr,"r=%g T=%g\n",c,radial_temp[i]);
    }*/
				       
    tbl = E->viscosity.zlith;
//    tbl = 0.86;
    tbl = E->sphere.rcomp;
//    bbl = E->sphere.rcomp;
    bbl = E->sphere.ri + 0.0223;
//    bbl = E->sphere.ri + 0.14;
    inttemp = 0.99;
    adgrad = 0.0;
    adtemp = inttemp + adgrad*(E->sphere.ro - E->sphere.ri)/E->sphere.ro;

		    
    noy=E->mesh.noy;  
    noz=E->mesh.noz;  
    nox=E->mesh.nox;  

    para1 = E->data.surf_temp*E->data.ref_temperature + 0.4*E->data.ref_temperature; 

    tbase = (para1 - E->data.surf_temp*E->data.ref_temperature)/E->data.ref_temperature;
    tbase1 = (para1 + 200 - E->data.surf_temp*E->data.ref_temperature)/E->data.ref_temperature;

    tbase = 0.5;

        mm = E->convection.perturb_mm[0];
        ll = E->convection.perturb_ll[0];

/*
        noz2 = (noz-1)/2+1;
        con = (noz-1)/(E->sphere.ro-E->sphere.ri);
*/
        con = E->convection.perturb_mag[0];

    if (E->control.restart==0 )    { 

          for(i=1;i<=noy;i++)
            for(j=1;j<=nox;j++)
              for(k=1;k<=noz;k++)  {
                node=k+(j-1)*noz+(i-1)*nox*noz;
                t1=E->X[1][node];
                r1=E->X[2][node];

                E->T[node] = 0.0;
                E->C[node] = 0.0;

/*
                if (k==noz2)  {
                  E->T[node] = con*modified_plgndr_a(ll,mm,t1);
                  }
*/
                E->T[node] = E->sphere.ri/(-E->sphere.ro+E->sphere.ri)*
                   (1.0-1.0/r1) + con*modified_plgndr_a(ll,mm,t1)*
                   sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));
                
		/*
              E->T[node] = tbase
                    + con*modified_plgndr_a(ll,mm,t1)*
                   sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));
		   */

		/*

              if (r1>=E->viscosity.zlith)
                E->T[node] = (r1-E->viscosity.zlith)/(E->sphere.ro-E->viscosity.zlith)*(0.0-tbase)+ tbase
               + con*modified_plgndr_a(ll,mm,t1)*sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));
	      else  if (r1<=(E->sphere.ri+0.03))
                E->T[node] = (r1-E->sphere.ri)/(0.03)*(tbase-1.0)+ 1.0 +con*modified_plgndr_a(ll,mm,t1)*sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));
	      else 
                E->T[node] = tbase
                    + con*modified_plgndr_a(ll,mm,t1)*
                   sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));
		   */

		if (r1>=tbl) {
		  T0 = E->TB[2][((node-1)/E->mesh.noz)*E->mesh.noz];
                  E->T[node] = T0 + (E->sphere.ro-r1)/(E->sphere.ro-tbl) * (inttemp-T0) + con*modified_plgndr_a(ll,mm,t1)*sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));
		}

                else if (r1 > bbl)
                  E->T[node] = inttemp + adgrad*(E->sphere.ro-r1)/E->sphere.ro + con*modified_plgndr_a(ll,mm,t1)*sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));

                else
                  E->T[node] = adtemp + (1.0-adtemp) * (bbl-r1)/(bbl-E->sphere.ri) + con*modified_plgndr_a(ll,mm,t1)*sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));

//		E->T[node] = radial_temp[k] + con*modified_plgndr_a(ll,mm,t1)*sin(M_PI*((r1-E->sphere.ri)/(E->sphere.ro-E->sphere.ri)));
		
		
              E->node[node] = E->node[node] | (INTX | INTZ | INTY);

              }   /* close the loop for node */

      if (E->control.composition) {
         if (!(strcmp(E->control.comp_adv_method,"field")==0))
            convection_initial_markers(E);
         else if ((strcmp(E->control.comp_adv_method,"field")==0))  {
            for(node=1;node<=E->mesh.nno;node++)  {
                  t1=E->X[1][node];
                  r1=E->X[2][node];
                  if (r1<=E->sphere.rcomp)
                     E->C[node] = 1.0;
                  else
                     E->C[node] = 0.0;
                  }
            }
         }
			
       }
	
   else if (E->control.restart==1)  {
      process_restart_tc(E,E->mesh.levmax);

      if (E->control.composition && !(strcmp(E->control.comp_adv_method,"field")==0)) {
          process_restart_mk(E); 
      }

     }



   temperatures_conform_bcs(E);


    thermal_buoyancy(E);


    return; 
    }

//  ------------ setup initial traces. could be done differently, for example
// using random positions.

 void convection_initial_markers(E)
   struct All_variables *E;
  {
    int el,i,j,k,p,node,ii,jj;
    float half_dist, init_height,dx,dr;
    void get_C_from_markers();
    double drand48();

    E->advection.element[0] =(int *)malloc((E->mesh.nel+1)*sizeof(int));
    E->advection.element[1] =(int *)malloc((E->mesh.nel+1)*sizeof(int));
    E->advection.element_prev[0] =(int *)malloc((E->mesh.nel+1)*sizeof(int));
    E->advection.element_prev[1] =(int *)malloc((E->mesh.nel+1)*sizeof(int));
	
    node = 0;
    p = pow((double)E->advection.markers_per_ele,(double)(1.0/E->mesh.dof));
    for (el=1;el<=E->mesh.nel;el++)  {
      dx = (E->X[1][E->ien[el].node[3]] - E->X[1][E->ien[el].node[1]])/p;
      dr = (E->X[2][E->ien[el].node[3]] - E->X[2][E->ien[el].node[1]])/p;
      for (i=1;i<=p;i++)
      for (j=1;j<=p;j++)  {
        node ++;
        E->XMC[1][node] = E->X[1][E->ien[el].node[1]] + dx*(i-0.5);
        E->XMC[2][node] = E->X[2][E->ien[el].node[1]] + dr*(j-0.5);
        E->CElement[node] = el;
        if (E->XMC[2][node]>E->sphere.rcomp){
              E->C12[node] = 0;
              E->C12f[node] = 0.0;
              }
        else {
              E->C12[node] = 1;
              E->C12f[node] = 1.0;
              }
        }
      }

    E->advection.markers = node;

    get_C_from_markers(E,E->C,E->CElement);
 return;
  }

/* ====================================================================== */

void process_restart_tc(E,lev)
    struct All_variables *E;
   int lev;
{
  int fileid[20];
  int i,j,k,ii,size2;
  char input_s[200],output_file[255],in_file[255];
  FILE *fp;
  float t1,t2;
  float t0,t3,t4,t5,t6,t7;

  sprintf(output_file,"%s/temp_comp.%d",E->convection.old_T_file,E->monitor.solution_cycles);

  fp = fopen(output_file,"r");

  fgets(input_s,200,fp);
  //sscanf(input_s,"%d %d %g",&i,&E->advection.timesteps,&E->monitor.elapsed_time);
  sscanf(input_s,"%d %d %g %g",&i,&E->advection.timesteps,&E->monitor.elapsed_time,&E->monitor.deltah);
  if (E->control.composition)   {
    for (i=1;i<=E->mesh.NNO[lev];i++)   {
      fgets(input_s,200,fp);
      sscanf(input_s,"%g %g %g %g",&E->T[i],&E->C[i],&t1,&t2);
      E->U[E->id[i].doff[1]]=t1;
      E->U[E->id[i].doff[2]]=t2;
      }
    for (i=1;i<=E->mesh.NEL[lev];i++)   {
      fgets(input_s,200,fp);
      sscanf(input_s,"%g",&t1);
      E->P[i] = t1;
      }
  }
  else {
    for (i=1;i<=E->mesh.NNO[lev];i++)   {
      fgets(input_s,200,fp);
      sscanf(input_s,"%g %g %g",&E->T[i],&t1,&t2);
      E->U[E->id[i].doff[1]]=t1;
      E->U[E->id[i].doff[2]]=t2;
      E->C[i]=0;
      }
    for (i=1;i<=E->mesh.NEL[lev];i++)   {
      fgets(input_s,200,fp);
      sscanf(input_s,"%g",&t1);
      E->P[i] = t1;
      }
    }

  E->advection.timesteps = E->monitor.solution_cycles;

  fclose(fp);

	/* Read in bulk composition and melt total */
  sprintf(output_file,"%s/ave.%d",E->convection.old_T_file,
					E->monitor.solution_cycles);

	fp = fopen(output_file,"r");

  fgets(input_s,200,fp);
  fprintf(stdout,"Reading from %s\n",input_s);
	sscanf(input_s,"%d %d %g %g %g %g %g %g %g %g %g\n",&t0,&t1,&t2,&t3,&t4,&t4,&t5,&t6,&t7,&E->Total.melt_prod,&E->Total.bulk_comp);

	fclose(fp);

	/* Read in slice.melt */
  sprintf(output_file,"%s/esurf.%d",E->convection.old_T_file,
					E->monitor.solution_cycles);

	fp = fopen(output_file,"r");

  for (i=1;i<=E->mesh.elx;i++)  {
     fgets(input_s,200,fp);
     sscanf(input_s,"%g %g %g",&t0,&E->slice.melt[i],&t1);
     }

  fclose(fp);

  return;
  }

/********************************************************
 * function read_tidal_heating                 		*
 *							*
 * Parameters: none					*
 *   	Reads in a table of tidal heating values on a 	*
 *	regular grid, output by e.g. Mathematica	*
 *							*
 * Returns: void					*
 *   	Gets the tidal heating for each Citcom element. *
 *   							*
 ********************************************************/

void read_tidal_heating(E)
    struct All_variables *E;
{
  int i,j,k;		/* Nodal counters */
  int ll,ul,ur,lr; 	/* Lower and upper left and right nodes for interp. */
  int numtheta,numr;	/* Number of regular gridpoints for tidal file */
  int e,ee,ends;					/* Element counter */
  char input_s[200],output_file[255],in_file[255];
  FILE *fp,*fp1,*fp11;
  float *theta,*r,*tidal_reg; 		/* Data from tidal heating input file */
  float *tidal_node, *tidal_elem;	/* Tidal Heating on Citcom grid */
  float A_ll, A_ul, A_lr, A_ur, A_node, A_elem;		/* Nodal "areas" */
  float temp;
  float totaltidal,totalvol;	/* Total tidal heating in mantle */
  float theta_n, a, c_h, d;	/* Shear heating distribution */

  ends = enodes[E->mesh.nsd];
	
  numtheta=91;
  numr=21;

  int input_int(char*, int*, char*);
  input_int("numtheta",&numtheta,"91");
  input_int("numr",&numr,"21");
  theta  = (float *) malloc((numtheta*numr+1)*sizeof(float));
  r  = (float *) malloc((numtheta*numr+1)*sizeof(float));
  tidal_reg = (float *) malloc((numtheta*numr+1)*sizeof(float));
  tidal_node = (float *) malloc((E->mesh.nno+1)*sizeof(float));
  tidal_elem = (float *) malloc((E->mesh.nel+1)*sizeof(float));

  /* Read in Tidal Heating on regular grid */

  sprintf(in_file,"%s/h2D.dat",E->control.data_file);
  fprintf(stderr,"Reading in tidal heating from '%s'\n",in_file);
  fp = fopen(in_file,"r");

  for (i=1;i<=(numtheta*numr);i++)   {
    fgets(input_s,200,fp);
    sscanf(input_s,"%g %g %g",&theta[i],&r[i],&tidal_reg[i]);
    /* Nondimensionalize inputs */
    tidal_reg[i] *= E->sphere.ro_dim*E->sphere.ro_dim
								/(E->data.density*E->data.Cp*E->data.DeltaT*E->data.therm_diff);
    }

  fclose(fp);
  fflush(fp);

  /* Interpolate from regular grid onto Citcom elements */

  for (e=1;e<=E->mesh.nel;e++) {  /* For each Citcom node */

 	/*
	 * Scan the regular grid to find nearest neighbors
	 * Probably there is a more efficient way to do this,
	 * but it's only done once, so who cares.
	 */
//	fprintf(stderr,"e %d %g %g\n",e,E->eco[e].centre[1],E->eco[e].centre[2]);
    for (i=1; i<=(numtheta*numr); i++) {
      /* If the citcom elem is between two adjacent regular nodes horiz.*/
      if((theta[i] <= E->eco[e].centre[1]) && (theta[i+numr] >= E->eco[e].centre[1])) {
//	fprintf(stderr,"Theta %d %g %g %g\n",i,theta[i],r[i],E->eco[e].centre[2]);
        /* and vertically */
        if((r[i] <= E->eco[e].centre[2]) && (r[i+1] >= E->eco[e].centre[2])) {
//	  fprintf(stderr,"r %d %g\n",i,r[i]);
	  ll = i;
	  ul = i+1;
	  lr = i+numr;
	  ur = i+numr+1;
 	  }
        }
      }
// fprintf(stderr,"ok4 %d\n",e);

    /* 
     * Average the heating values at the four regular elements adjacent to 
     * the Citcom element in question, weighted by their proximity.
     */
//fprintf(stderr,"%d %d %d %d\n",ll,ul,lr,ur);
if (ll==17) {
fprintf(stdout,"%d %g %g\n",ll,theta[ll],r[ll]);
fprintf(stdout,"%d %g %g\n",ul,theta[ul],r[ul]);
fprintf(stdout,"%d %g %g\n",lr,theta[lr],r[lr]);
fprintf(stdout,"%d %g %g\n",ur,theta[ur],r[ur]);
}

    /* Areas of quadrangle formed by regular and Citcom elements */
    A_ll = (theta[ll] - E->eco[e].centre[1]) * (r[ll] - E->eco[e].centre[2]);
    A_ul = (theta[ul] - E->eco[e].centre[1]) * -(r[ul] - E->eco[e].centre[2]);
    A_lr = -(theta[lr] - E->eco[e].centre[1]) * (r[lr] - E->eco[e].centre[2]);
    A_ur = -(theta[ur] - E->eco[e].centre[1]) * -(r[ur] - E->eco[e].centre[2]);

    /*Total Area of quadrangle formed by adjacent elements*/
    A_elem = A_ll + A_ul + A_lr + A_ur;

    /* Weighting is proportional to the area of the opposite regular node */
    /* Larger area means node is farther away.  Lever rule. */
    tidal_elem[e] = (tidal_reg[ll]*A_ur + tidal_reg[lr]*A_ul 
			+ tidal_reg[ul]*A_lr + tidal_reg[ur]*A_ll) / A_elem;

if (ll==17) {
fprintf(stdout,"%d %g %g\n",ll,tidal_reg[ll],A_ll);
fprintf(stdout,"%d %g %g\n",ul,tidal_reg[ul],A_ul);
fprintf(stdout,"%d %g %g\n",lr,tidal_reg[lr],A_lr);
fprintf(stdout,"%d %g %g\n",ur,tidal_reg[ur],A_ur);
fprintf(stdout,"%d %g %g %g\n\n",e,tidal_elem[e],E->eco[e].centre[1],E->eco[e].centre[2]);
}
    }

  /* Project from nodes to elements */
  totaltidal = 0.0;
  totalvol = 0.0;

  for (e=1;e<=E->mesh.nel;e++)  {

    E->heating_tidal[e] = tidal_elem[e];

    totaltidal += E->heating_tidal[e]*E->eco[e].area;
    totalvol += E->eco[e].area;

    }

  /* Add Shear heating near surface */
  if(E->control.shear_heating) {
    /* Disk-shaped region centered on S pole */
    theta_n = 145.0 * (M_PI / 180.0);  	/* Northernmost extent */
    a = 0.28;				/* Broadness of gaussian */
    d = 5.0e3;				/* Thickness of heated layer in m*/
    c_h = 0.273 / d;			/* Heating rate at pole (W/m^3) */

    /* Nondimensionalize c_h */
    c_h *= E->sphere.ro_dim*E->sphere.ro_dim
						/ (E->data.density*E->data.Cp*E->data.DeltaT*E->data.therm_diff);

    for (e=1;e<=E->mesh.nel;e++) {
      ee = e + (e-1) / (E->mesh.elz); 
      if ((E->X[1][ee] >= theta_n) 
			&& (E->X[2][ee] >= (1.0 - d/E->sphere.ro_dim))) {
	/* Gaussian distribution */
	/*E->heating_shear[e] = (2.0/(a*a)) 
		* exp(-(M_PI-E->X[1][ee])*(M_PI-E->X[1][ee]) / a*a);
	*/
	/* Sinusoidal Distribution */
	E->heating_shear[e] = c_h 
		* cos(M_PI/2.0 * (M_PI-E->X[1][ee]) / (M_PI-theta_n)) 
		* cos(M_PI/2.0 * (M_PI-E->X[1][ee]) / (M_PI-theta_n));
      }
      else E->heating_shear[e] = 0.0;
    }
 
  }

  /* Output tidal and shear heating to file */

  sprintf(output_file,"%s/heating_t.%d",E->control.data_file,E->monitor.solution_cycles);
  fp11 = fopen(output_file,"w");
	fprintf(stdout,"%f %f\n",totaltidal,totalvol);
  for (e=1;e<=E->mesh.nel;e++){
    fprintf(fp11,"%g %g %g %g\n",E->eco[e].centre[1],E->eco[e].centre[2],
	E->heating_tidal[e],E->heating_shear[e]);
    }
fflush(fp11);

/*fprintf(stderr,"ok6 \n");*/

    free((void *) theta );  
    free((void *) r );  
    free((void *) tidal_reg );  
    free((void *) tidal_node );  
    fclose(fp11);
/*fprintf(stderr,"ok7 \n");*/


  return;
}


/********************************************************
 * function read_impacts			                 					*
 *																											*
 * Parameters: 																					*
 *	num_impacts	-- The total number of impacts					*
 *																											*
 * Returns: void																				*
 *   	Reads in the parameters from each impact from a	 	*
 *		file specified in the input_file									*
 *   																										*
 ********************************************************/

void read_impacts(E,num_impacts)
	struct All_variables *E;
	int num_impacts;
{
	char impact_file[255], input_s[1000];
	const float rad=180.0/M_PI; /* radian to degree convesion */
	float start;
	float ksi,C,S;							/* Shock parameters */
	int i; 											/* counter */

	FILE *fp;

  /* Allocate transient depth pointer */
  E->impacts.H_t     = (float *)malloc((E->mesh.nsf+2)*sizeof(float));

	fprintf(stdout,"Read Impacts: %d\n",num_impacts);
	fprintf(E->fp,"Read Impacts\n");

  int input_int(char*, int*, char*);
  int input_string();
  input_int("impact_heating_from_file",&(E->impacts.heat_from_file),"0");
  if(E->impacts.heat_from_file)
      input_string("impact_heating_file",E->impacts.heating_file);
	fprintf(stdout,"impact heating file %s\n",E->impacts.heating_file);

  /* Open file containing impact parameters*/
	sprintf(impact_file,"%s/impact_file",E->control.data_file);
	fprintf(stdout,"%s\n",impact_file);
  fprintf(E->fp,"%s\n",impact_file);
	fp=fopen(impact_file,"r");

  fprintf(E->fp,"Starting age %f My\n",E->control.restart_age);

	/* Determine temperature increase as function of impact velocity */
  int input_float(char*, float*, char*);
  input_float("impact_velocity",&(E->impacts.v),"15.0e3");

	/* From Watters et al. (2008) */
	S = 1.25;
	C = 7.4e3;
	
	ksi = C/(S*0.5*E->impacts.v);
	E->impacts.dT = (E->impacts.v*E->impacts.v/4.0)/(2.0*E->data.Cp)
        							* (1.0 - 2.0*(ksi - ksi*ksi*log(1.0 + 1.0/ksi)));


	/* Read in the impact data */
  fgets(input_s,1000,fp);  /* Header Line */

  for(i=0;i<num_impacts;i++) {
    fgets(input_s,1000,fp);
    sscanf(input_s,"%f %f %f %f",&(E->impacts.t[i]),&(E->impacts.f[i]),
				&(E->impacts.th[i]),&(E->impacts.size[i]));

    fprintf(E->fp,"%f %f %f %f %f\n",E->impacts.t[i],E->impacts.th[i],
			E->impacts.f[i],E->impacts.v,E->impacts.size[i]);
    fprintf(stderr,"%f %f %f %f %f\n",E->impacts.t[i],E->impacts.th[i],
			E->impacts.f[i],E->impacts.v,E->impacts.size[i]);
  }

  fclose(fp);

	/*Nondimensionalize impact data*/
  E->impacts.dT *= 1.0/E->data.DeltaT;
  start = E->control.restart_age/1.0e3;

  for(i=0;i<E->impacts.number;i++) {
    E->impacts.t[i] = (start - E->impacts.t[i]) * 1.0e9;
    E->impacts.t[i] *= (86400.0*365.25) * E->data.therm_diff
                        / (E->sphere.ro_dim*E->sphere.ro_dim);
    E->impacts.t[i] += E->monitor.elapsed_time;
    E->impacts.th[i] = (90.0 - E->impacts.th[i]) / rad;
    E->impacts.f[i] *= 1.0/rad; 
    E->impacts.flag[i] = 0;
    
		/* If time of impact is before restart time, set size to 0 */
    if((E->impacts.t[i] < E->monitor.elapsed_time) || (E->impacts.t[i] < 0.0) ){
        E->impacts.size[i] = 0.0;
        E->impacts.flag[i] = 1;
        }

		fprintf(E->fp,"%e %e %f %f %f %d\n",start,E->impacts.t[i],E->impacts.th[i],E->impacts.f[i],E->impacts.size[i],E->impacts.flag[i]);
	}

}

/********************************************************
 * define_solidus                                       *
 *      function to define the solidus vs. depth in the *
 *      mantle.  Taken from Reese et al. (2002)         *
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 ********************************************************/

void define_solidus(struct All_variables *E)
{
  int k;        /* Radial node */
  float z;      /* Depth (m) */
  float P,Pc;   /* Pressure (GPa) */

  Pc = 10.7;    /* P at which solidus fn changes */

	fprintf(E->fp,"Solidus: \n");
	for(k=1;k<=E->mesh.noz;k++) {
    z = (1.0 - E->X[2][k])*E->sphere.ro_dim;
    P = E->data.density*E->data.grav_acc*z/1.0e9;
    /* Reese et al. 2002*/
    /*if (P <= Pc)
      E->solidus[k] = 1374.0 + 130.0*P - 5.6*P*P;
    else
      E->solidus[k] = 2017.0 + 10.0*P;
    */

    /* Katz et al. 2003 (scaled to Kelvins) */
    /*
    E->solidus[k] = 1358.7 + 132.9*P -5.1*P*P;
    E->lherzliq[k] = 1748.0 + 80.0*P -3.2*P*P;
    E->liquidus[k] = 2053.0 + 45.0*P -2.0*P*P;
    */

    /* Water ice, pure */
    E->solidus[k] = 273.0;  /* P-dep v. small on Enc.*/
    E->liquidus[k] = 273.0;  /* P-dep v. small on Enc.*/

    /* Nondimensionalize */
    E->solidus[k] = E->solidus[k]/E->data.DeltaT - E->data.surf_temp;
//    E->lherzliq[k] = E->lherzliq[k]/E->data.DeltaT - E->data.surf_temp;
    E->liquidus[k] = E->liquidus[k]/E->data.DeltaT - E->data.surf_temp;
      
	/*	fprintf(E->fp,"%d %e %f %f %f %f \n",k,z,P,E->solidus[k],
																		E->lherzliq[k],E->liquidus[k]);
                                    */
		fprintf(E->fp,"%d %e %f %f %f\n",k,z,P,E->solidus[k],E->liquidus[k]);
  }

}

/****************************************************************
 * despin																												*
 * 																															*
 * Computes the stresses and heating due to despinning of a 		*
 * body based on the orbital evolution and mechanical 					*
 * properties.																									*
 ****************************************************************/

void despin(E)
  struct All_variables *E;
  {
   
	double qpi;  
	complex long double c_rigidity;	/* Viscoelastic "rigidity" */
	complex double c_k2;					/* Love Number k2 */
  double mu = E->data.rigidity;
  double eta = E->data.ref_viscosity;
  double omega;
  long double mean_motion;		/* Initial orbital frequency */
  long double semimajor_init;	/* Initial orbital distance */
  long double semimajor_final;	/* Final orbital distance */
  double drotdt;		/* Despinning rate */
	double Ls1, Ls2, Lo1, Lo2;
	double Es1, Es2, Eo1, Eo2;
	double vol;


	qpi = 3.141592653589793;
  E->control.despun = 0;	/* Flag for completion of despinning */

	E->data.mass = ((E->data.density_core-E->data.density) * pow(E->sphere.ri,3) 
		  + E->data.density) * pow(E->sphere.ro_dim,3) * 4.0*(qpi)/3.0;
	vol = pow(E->sphere.ro_dim,3) * 4.0*(qpi)/3.0;

																													/* Dimensional */
  E->data.moi *= (E->data.mass*E->sphere.ro_dim*E->sphere.ro_dim); 

  fprintf(stderr,"PI %.15f, Volume %.15e\n",qpi,vol);
  fprintf(stderr,"densities %g, %g;  ro %g ri %g \n mass: %.15e kg  moi: %.15e\n",
		E->data.density_core,E->data.density,E->sphere.ro_dim,E->sphere.ri,
		E->data.mass,E->data.moi);

/* 
 * Compute total dissipated energy 				
 * We know the moment of inertia, initial and final rotation 	
 * rates, mass, and final orbital semimajor axis.  Really, 
 * these are time-dependent, and we should evolve them too.
 */

 /* First, solve for initial semimajor axis using conservation 
    of angular  momentum */

  semimajor_final = powl((E->data.grav_const * E->data.mass_primary)
														/(E->data.rot_final*E->data.rot_final),(1.0/3.0));

  Ls1 = E->data.moi*E->data.rot_init;
  Ls2 = E->data.moi*E->data.rot_final;
	Lo2 = E->data.mass*E->data.rot_final*semimajor_final*semimajor_final;
	
  Lo1 = Ls2 + Lo2 - Ls1;

   semimajor_init = powl((Lo1/E->data.mass),2) 
		 											/ (E->data.grav_const * E->data.mass_primary);
	

/*  semimajor_init = pow((E->data.moi*(E->data.rot_final - E->data.rot_init)
																															/ E->data.mass + 
					E->data.semimajor_axis*E->data.semimajor_axis*E->data.rot_final),2) 
	  								/ (E->data.grav_const * E->data.mass_primary);
*/
	E->data.semimajor_axis = semimajor_init;

  /* Next, solve for initial non-synchronous mean motion using 
     Kepler's 3rd Law */

  E->data.mean_motion = sqrtl((E->data.grav_const * E->data.mass_primary) / 
								  powl(semimajor_init,3));

  /*
	fprintf(stderr,"Ls1 %lg Ls2 %g Lo1 %g Lo2 %g\n",Ls1,Ls2,Lo1,Lo2);
	fprintf(stderr,"Ls1 - Ls2 %g Lo1 - Lo2 %g\n",Ls1-Ls2,Lo1-Lo2);
  */

  /* Finally, find dissipated energy using conservation of energy */

	Es1 =  E->data.moi*E->data.rot_init*E->data.rot_init /2.0;
	Es2 =  E->data.moi*E->data.rot_final*E->data.rot_final /2.0;
	Eo1 = - E->data.grav_const*E->data.mass_primary*E->data.mass 
				/ (2.0*semimajor_init);
	Eo2 = - E->data.grav_const*E->data.mass_primary*E->data.mass 
				/ (2.0*semimajor_final);

  E->data.Ediss_total = ( E->data.moi*(E->data.rot_init*E->data.rot_init 
								  					- E->data.rot_final*E->data.rot_final) 
	  										- E->data.grav_const*E->data.mass_primary*E->data.mass 
		  						* (1.0/semimajor_init - 1.0/semimajor_final) )/2.0;

	E->data.Ediss = 0.0; /* Initialize running total */

	fprintf(stderr,"rot2 %g rot1 %g M %g \n af %.10Lg ai %.7Lg n %Lg E %g\n",
		E->data.rot_final,E->data.rot_init,E->data.mass_primary,
		semimajor_final,semimajor_init,E->data.mean_motion,E->data.Ediss_total);
	fprintf(stderr,"Es1 %g Es2 %g Eo1 %g Eo2 %g\n",Es1,Es2,Eo1,Eo2);
	fprintf(stderr,"Es1 - Es2 %g Eo1 - Eo2 %g\n",Es1-Es2,Eo1-Eo2);

	fprintf(E->fp,"rot2 %g rot1 %g M %g \n af %.10Lg ai %.7Lg n %Lg E %g\n",
		E->data.rot_final,E->data.rot_init,E->data.mass_primary,
		semimajor_final,semimajor_init,E->data.mean_motion,E->data.Ediss_total);
	fprintf(E->fp,"Es1 %g Es2 %g Eo1 %g Eo2 %g\n",Es1,Es2,Eo1,Eo2);
	fprintf(E->fp,"Es1 - Es2 %g Eo1 - Eo2 %g\n",Es1-Es2,Eo1-Eo2);

/* 
 * Compute initial despinning rate and despinning timescale assuming 
 * constant rate.
 */

  omega = E->data.rot_init;
	fprintf(stderr,"omega %.8g mu %.8g eta %.8g\n",omega,mu,eta);
  c_rigidity = omega*eta*mu/(mu*mu + omega*omega*eta*eta);
  c_rigidity *= (omega*eta + I*mu);
  c_k2 = 1.5 / (1.0 + (19.0) * c_rigidity 
											/ (2*E->data.grav_acc*E->data.density*E->sphere.ro_dim));

  drotdt = (3.0 * -cimag(c_k2) * E->data.grav_const * E->data.mass_primary 
							* E->data.mass_primary * pow((E->sphere.ro_dim),5)) 
							/ (E->data.moi * pow(E->data.semimajor_axis,6));

	fprintf(stderr,"rigidity: %.7e + %.7e i Pa, k2: %.7e + %.7e i, \nrate: %.8g /s\n",creal(c_rigidity),cimag(c_rigidity),creal(c_k2),cimag(c_k2),drotdt);
  
  E->data.rotation = E->data.rot_init;
  E->data.despin_rate = drotdt;

  /* Initially Hydrostatic planet */
	E->data.flattening = 1.25*E->data.rot_init*E->data.rot_init*E->sphere.ro_dim 
													/ E->data.grav_acc;

	/* Assuming constant despin rate ... */

  E->control.despin_timescale = omega / drotdt;
  E->power_despin_global = (E->data.Ediss_total / E->control.despin_timescale);
  
  /* Volumetric Average */
  E->power_despin_ave = E->power_despin_global / ((4.0*M_PI/3.0)
                    *pow(E->sphere.ro_dim,3) * (1.0 - pow(E->sphere.ri,3)) ); 
  fprintf(stderr,"Dissipated energy: %.7g J %.7g W %.7g W/m^3\n",
							E->data.Ediss_total, E->power_despin_global, E->power_despin_ave);
  fprintf(stderr,"Despin Timescale: %.7g s = %.7g y\n",
							E->control.despin_timescale,E->control.despin_timescale/3.15e7);

  /* Nondimensionalize */
  E->power_despin_global *= 1.0/(E->data.density*E->data.Cp*E->data.DeltaT
					*E->data.therm_diff*E->sphere.ro_dim);
  E->power_despin_ave *= (3.0/(4.0*M_PI*(1.0 - pow(E->sphere.ri,3))))
		* E->power_despin_global;
  E->control.despin_timescale *= E->data.therm_diff 
				/ (E->sphere.ro_dim*E->sphere.ro_dim);
  } 



/* ====================================================================== */

void process_restart_mk(E)
  struct All_variables *E;
{
 int fileid[20];
 int i,j,k,ii,size2;
 int t0,t1;
 char input_s[200],output_file[255],in_file[255];
 FILE *fp;
 void p_to_centres();
 void get_C_from_markers();

    E->advection.element[0] =(int *)malloc((E->mesh.nel+1)*sizeof(int));
    E->advection.element[1] =(int *)malloc((E->mesh.nel+1)*sizeof(int));
    E->advection.element_prev[0] =(int *)malloc((E->mesh.nel+1)*sizeof(int));
    E->advection.element_prev[1] =(int *)malloc((E->mesh.nel+1)*sizeof(int));
 sprintf(output_file,"%s/traces.%d",E->convection.old_T_file,E->monitor.solution_cycles);

 fp = fopen(output_file,"r");

   fgets(input_s,200,fp);
   sscanf(input_s,"%d %d %g",&E->advection.markers,
	 				&E->advection.timesteps,&E->monitor.elapsed_time);

   for (i=1;i<=E->advection.markers;i++)  {
     fgets(input_s,200,fp);
     /*sscanf(input_s,"%g %g %d %d",&E->XMC[1][i],&E->XMC[2][i],
						&E->CElement[i],&E->C12[i]);*/
     sscanf(input_s,"%g %g %d %g",&E->XMC[1][i],&E->XMC[2][i],
						&E->CElement[i],&E->C12f[i]);
     }

   for (i=1;i<=E->mesh.nel;i++)  {
     fgets(input_s,200,fp);
     sscanf(input_s,"%g %d %d",&E->CE[i],&t0,&t1);
     }

  E->advection.timesteps = E->monitor.solution_cycles;

	p_to_centres(E,E->C,E->CE,E->mesh.levmax);
  get_C_from_markers(E,E->C,E->CElement);
     //   	reallocate_markers(E);
  	//			get_markers_from_C(E,E->CE,E->CElement);

  fclose(fp);

  return;
  }


/********************************************************
 * adjust_model_domain                                  *
 *      function to adjust the model domain, including  *
 *      grids, global shape function, and interpolation *
 *      of relevant quantities onto the new grid in the *
 *      case where the ice shell thickens or thins      *
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *                                                      *
 * Returns                                              *
 *      none                                            *
 ********************************************************/
void adjust_model_domain(E)
  struct All_variables *E;
{
  static int visits=0;
  
  void allocate_new_vars();
  void update_node_locations();
  void interpolate_vars();
  void mass_matrix();

  fprintf(stderr,"Adjusting Model Domain, deltah = %f m\n",
        E->sphere.deltarb*E->sphere.ro_dim);
  fprintf(E->fp,"Adjusting Model Domain, deltah = %f m\n",
        E->sphere.deltarb*E->sphere.ro_dim);

  if(visits==0) {
    allocate_new_vars(E);  /* Allocates new P, T, visc arrays. */
    visits++;
    fprintf(stderr,"Allocate New Vars\n");
  }
  update_node_locations(E);   /* Defines new radial coords for nodes */
  fprintf(stderr,"Update Node Locations\n");
  mass_matrix(E);             /* Rebuild the mass matrix and compute new */
                              /*   element coordinates */
  fprintf(stderr,"Rebuild Mass Matrix\n");
  interpolate_vars(E);    /* Interpolate T, C, onto new grid */
  fprintf(stderr,"Interpolate Vars\n");

  E->sphere.deltarb = 0.0;  /* Reset deltarb */

}


/********************************************************
 * allocate_new_vars                                    *
 *      function to allocate memory for arrays for      *
 *      pressure, temperature, and viscosity            *
 *      interpolated from the old grid onto the new     *
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *                                                      *
 * Returns                                              *
 *      none                                            *
 ********************************************************/
void allocate_new_vars(E)
  struct All_variables *E;
{
   
  int i,j;           /* Counters */

  /* Allocate memory for key arrays on new grid */   
  for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)   
    E->newvars.XX[i] = (float *)  malloc((E->mesh.NNO[i]+1)*sizeof(float));
    
  E->newvars.XP = (float *) malloc ((E->mesh.noz+1)*sizeof(float));
  E->newvars.oldX = (float *) malloc ((E->mesh.nno+1)*sizeof(float));
  E->newvars.oldecocentre = (float *) malloc ((E->mesh.nel+1)*sizeof(float));

  E->newvars.C = (float *) malloc ((E->mesh.nno+1)*sizeof(float));
  E->newvars.Fm = (float *) malloc ((E->mesh.nno+1)*sizeof(float));
  E->newvars.T = (float *) malloc ((E->mesh.nno+1)*sizeof(float));
  E->newvars.P = (double *) malloc ((E->mesh.npno+1)*sizeof(double));
  E->newvars.heating_tidal = (float *) malloc ((E->mesh.nel+1)*sizeof(float));
  E->newvars.tidal_visc = (float *) malloc ((E->mesh.nel+1)*sizeof(float));
  E->newvars.mat = (int *) malloc((E->mesh.nel+2)*sizeof(int));
  
  /* Ones originally allocated in allocate_velocity_vars */
  E->newvars.U = (double *) malloc((E->mesh.nsd*E->mesh.nnov+1)*sizeof(double));

    /* Initialize */
    for(i=1;i<=E->mesh.nno;i++)
      E->newvars.T[i] = 0.0;

    for(i=1;i<E->mesh.nel;i++)   {
      E->newvars.mat[i]=1;
      E->newvars.heating_tidal[i] = 0.0;
    }
  
    for(i=1;i<=E->mesh.npno;i++)
      E->newvars.P[i] = 0.0;
}


/********************************************************
 * interpolate_vars                                     *
 *      function to interpolate key variables           *
 *      temperature and compostion onto new grid.       *
 *      buoyancy is recomputed from new T, C.           *
 *                                                      *
 *      For now, pressure and velocity are zeroed out   *
 *      and recomputed from scratch in next call to     *
 *      general_stokes_solver as for a new case. We     *
 *      plan to replace with an interpolation later.    *
 *      This should not be a problem for Newtonian      *
 *      viscosity, but nonNewt. visc depends on prev.   *
 *      V, so we'll need interp for that one.           *
 *                                                      *
 *      Tidal heating is not being interpolated. This   *
 *      needs to be completely recomputed from new      *
 *      viscosity structure and must wait for next call *
 *      to TiRADE.                                      *
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *                                                      *
 * Returns                                              *
 *      none                                            *
 ********************************************************/
void interpolate_vars(E)
  struct All_variables *E;
{
  float rnew,rold,rold2;           /* radial pos. in old and new grid */    

  int i,i2,j,k;              /* Counters */
  int nox,noy,noz;         /* Nodes in each dir */
  int node,node2;          /* Node indices */

  void thermal_buoyancy();

  noy=E->mesh.noy;  
  noz=E->mesh.noz;  
  nox=E->mesh.nox;  

  /* By node */
  for(i=1;i<=noy;i++)
    for(j=1;j<=nox;j++) {
      for(k=2;k<=(noz-1);k++)  {  /* Don't interp boundary nodes */
        node=k+(j-1)*noz+(i-1)*nox*noz;
        rnew=E->X[2][node]; /* Radius of node in new grid */
        /* Find nearest nodes in old grid */
        rold=E->newvars.oldX[node];           /* Radius of node in old grid */

        if(E->sphere.deltarb > 0 )    /* If shell thickens */
          node2 = node-1;
        else                            /* If shell thins */
          node2 = node+1;
        rold2=E->newvars.oldX[node2];        /* Radius of nearest other *
                                              * node in old grid        */
                                          
        /* Temperature and Composition */
        E->newvars.T[node] = ( E->T[node]*(rnew - rold2) 
                             + E->T[node2]*(rnew - rold) ) / (rold - rold2);
        
        E->newvars.C[node] = ( E->C[node]*(rnew - rold2) 
                             + E->C[node2]*(rnew - rold) ) / (rold - rold2);

      }

      /* Bottom (k=1) must keep same BC */
      node = 1+(j-1)*noz+(i-1)*nox*noz;
      E->newvars.T[node] = E->T[node]; 
      E->newvars.C[node] = E->C[node]; 
      /* Surface (k=noz) does not change */
      node = noz+(j-1)*noz+(i-1)*nox*noz;
      E->newvars.T[node] = E->T[node]; 
      E->newvars.C[node] = E->C[node]; 

    } /* End of node loop */

  /* Copy new T and C into original arrays */  
  for(i=1;i<=E->mesh.nno;i++) {
    E->T[node] = E->newvars.T[node];
    E->C[node] = E->newvars.C[node];
  }

  /* Get buoyancy force (function of T and C) */
  thermal_buoyancy(E);

  /* Pressure */
  for(i=1;i<=E->mesh.npno;i++) {
    /* Just zero it out for now, like initial state. Ultimately, we will
       want to fix this. Should be able to interp as with T, but how
       to handle P at bottom boundary? Probably just go with P at nearest node a       above */
    E->P[i] = 0.0;

    //  rnew=E->eco[i].centre[2]; /* Radius of pnode in new grid */
    /* Find nearest pnodes in old grid */
    // rold=E->newvars.oldecocentre[i]; 

    // if(E->sphere.deltarb > 0 )    /* If shell thickens */
    //   i2 = i-1;
    //else                            /* If shell thins */
    //   i2 = i+1;
    //rold2=E->newvars.oldecocentre[i2]; 
  }

  /* Velocity */
    /* For now, we zero it out like P. We need to interp both U and P 
       or neither */
  for(i=0;i<E->mesh.neq;i++)
    E->U[i]=0.0;
  
  for(i=1;i<=E->mesh.nnov;i++)   {
	  E->V[1][i]=0.0;
	  E->V[2][i]=0.0;
	}

  if(E->mesh.dof==3)   {
 	  for(i=1;i<=E->mesh.nnov;i++)
	    E->V[3][i]=0.0;
  }
 
}
