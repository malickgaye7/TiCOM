/*****************************************
 *   CC  III  TTTTT   CC   OO   MM MM    *
 *  C     I     T    C    O  O  M M M    *
 *  C     I     T    C    O  O  M   M    *
 *   CC  III    T     CC   OO   M   M    *
 *                                       *  
 * Developed at CIT for COnvection in    *
 * the Mantle by Louis Moresi 1992-today *
 *                                       *
 * You are free to use this code but it  * 
 * is distrubuted as BeWare i.e. it does *
 * not carry any guarantees or warranties *
 * of reliability.                       *
 *                                       *
 * Please respect all the time and work  *
 * that went into the development of the *
 * code.                                 *  
 *                                       *
 * LM                                    *
 *****************************************/
/*   Functions which solve the heat transport equations using Petrov-Galerkin
     streamline-upwind methods. The process is basically as described in Alex
     Brooks PhD thesis (Caltech) which refers back to Hughes, Liu and Brooks.  */

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif
//#include <malloc.h>
//#include <sys/types.h>

#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif

#ifndef __ELEMENT_DEFINITIONS_H__
#define __ELEMENT_DEFINITIONS_H__
#include "element_definitions.h"
#endif

#ifndef __GLOBAL_DEFS_H__
#define __GLOBAL_DEFS_H__
#include "global_defs.h"
#endif

#ifndef COMPLEX_H
#define COMPLEX_H
#include <complex.h>
#endif

#ifndef STRING_H
#define STRING_H
#include <string.h>
#endif

#include "Parsing.h"

extern int Emergency_stop;

struct el { double gpt[9]; };

/* ============================================
   Generic adv-diffusion for temperature field.
   ============================================ */


void advection_diffusion_parameters(E)
     struct All_variables *E;

{

    void std_timestep();
    /* Set intial values, defaults & read parameters*/
    
    E->advection.temp_iterations = 2; /* petrov-galerkin iterations: minimum value. */
    E->advection.total_timesteps = 1; 
    E->advection.sub_iterations = 1;
    E->advection.last_sub_iterations = 1;
    E->advection.gamma = 0.5;
    E->advection.dt_reduced = 1.0;         

    E->monitor.T_maxvaried = 1.02;
 
    input_boolean("ADV",&(E->advection.ADVECTION),"on");
    E->advection.ADVECTION=1;
   
    input_int("visc_heating",&(E->control.visc_heating),"1");
    input_int("adi_heating",&(E->control.adi_heating),"1");
    input_int("latent_heating",&(E->control.latent_heating),"1");

    input_int("minstep",&(E->advection.min_timesteps),"1");
    input_int("maxstep",&(E->advection.max_timesteps),"1000");
    input_int("maxtotstep",&(E->advection.max_total_timesteps),"1000000");
    input_float("finetunedt",&(E->advection.fine_tune_dt),"0.9");
    input_float("fixed_timestep",&(E->advection.fixed_timestep),"0.0");
    input_int("adv_sub_iterations",&(E->advection.temp_iterations),"2,2,nomax");
    input_float("maxadvtime",&(E->advection.max_dimensionless_time),"10.0");
   
    input_float("sub_tolerance",&(E->advection.vel_substep_aggression),"0.005");  
    input_int("maxsub",&(E->advection.max_substeps),"25");
  
    input_float("liddefvel",&(E->advection.lid_defining_velocity),"0.01");
    input_float("sublayerfrac",&(E->advection.sub_layer_sample_level),"0.5");            
        input_float("Q0_enriched",&(E->control.Q0ER),"0.0");
    input_int("markers_per_ele",&(E->advection.markers_per_ele),"0");
	      
  /* allocate memory */

  return;
}

void advection_diffusion_allocate_memory(E)
     struct All_variables *E;

{ int i;

  E->Tdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
  for(i=1;i<=E->mesh.nno;i++) 
    E->Tdot[i]=0.0;

  if (E->control.composition)   {
    E->Cdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    for(i=1;i<=E->mesh.nno;i++)
      E->Cdot[i]=0.0;
    if (!(strcmp(E->control.comp_adv_method,"field")==0)) {
      E->advection.markers = E->advection.markers_per_ele*E->mesh.nel;
      for(i=1;i<=E->mesh.nsd;i++)   {
        E->VO[i] = (float *) malloc ((E->advection.markers+1)*sizeof(float));
        E->XMC[i] = (float *) malloc ((E->advection.markers+1)*sizeof(float));
        E->XMCpred[i] = (float *) malloc ((E->advection.markers+1)*sizeof(float));
        }
      E->C12 = (int *) malloc ((E->advection.markers+1)*sizeof(int));
      E->C12f = (float *) malloc ((E->advection.markers+1)*sizeof(float));
      E->CElement = (int *) malloc ((E->advection.markers+1)*sizeof(int));
      }

    }

return;
}

void PG_timestep_particle(E)
     struct All_variables *E;
{
    void timestep();
    void predictor();
    void corrector();
    void pg_solver();
    void remove_horiz_ave();
    void return_horiz_ave();
    void std_timestep();
    void temperatures_conform_bcs();
    void thermal_buoyancy();
    void Runge_Kutta();
    void Euler();
    void get_fixed_temp();
		void impact_heating();
		void p_to_centres();
		void reallocate_markers();
		void get_markers_from_C();

		float bulk_comp();
		float melting();
		float melting_by_node();
    float Tmax(),T_interior1;
		float F, C;
    int e,i,j,k,psc_pass,count,steps,iredo;
    int keep_going;

    float *T1, *Tdot1;
    static float *DTdot;
    FILE *fp;

		const int nel=E->mesh.nel;

    static int loops_since_new_eta = 0;
    static int been_here = 0;
    static int been_here0 = 0;
    static int on_off = 0;

   if (been_here++==0)  {
      DTdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
      }

    if (on_off==0)  {
      E->advection.timesteps++;
      std_timestep(E);
      E->advection.total_timesteps++;


      }


    if (on_off==1)  {

        Runge_Kutta(E,E->XMC,E->XMCpred,E->C,E->V,E->VO);
 			  C = bulk_comp(E,2);

				/* After Euler, R-L steps, bulk comp should be unchanged
				 * Numerically, this is tricky because tracers do not have
				 * volumes associated with them and elements are not the 
				 * same size.  This is a kludge to maintain same comp
				 * on a global scale.  Melting, of course, changes things.
				 */

				for(i=1;i<=E->mesh.nno;i++){
					if (C != 0)
						E->C[i] *= ( E->Total.bulk_comp/C);
				}
				
        for(e=1;e<=E->mesh.nel;e++){
					if (C != 0)
						E->CE[e] *= ( E->Total.bulk_comp/C);
				}

 			  C = bulk_comp(E,3);

			  //p_to_centres(E,E->C,E->CE,E->mesh.levmax);
			  E->Total.bulk_comp_prev = 0.0;
			  for(e=1;e<=E->mesh.nel;e++) {
			   	E->Total.bulk_comp_prev += E->CE[e] * E->eco[e].area;
				  E->CE_prev[e] = E->CE[e];
			  }
			
			  E->Total.bulk_comp_prev /= E->Total.vol;

    }

    else if (on_off==0)   {
			/*Store old values */	

			for(i=1;i<=E->mesh.nno;i++)
				E->C_prev[i] = E->C[i];

			for(e=1;e<=E->mesh.nel;e++) {
					E->CE_prev[e] = E->CE[e];
					E->advection.element_prev[0][e] = E->advection.element[0][e];
					E->advection.element_prev[1][e] = E->advection.element[1][e];
			}

		//	C = bulk_comp(E,10);
      Euler(E,E->XMC,E->XMCpred,E->C,E->V,E->VO);
			C = bulk_comp(E,1);

                                 /* update temperature    */
      predictor(E,E->T,E->Tdot,1);
      for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++)   {
          pg_solver(E,E->T,E->Tdot,DTdot,E->V,1.0,E->TB,E->node);
          corrector(E,E->T,E->Tdot,DTdot,1);
          }

      temperatures_conform_bcs(E,E->T);

      }

		/* Need to track melting at all times, not just impacts */
		/* Only do this during the final pass */

		if (on_off == 1) {
			if (E->control.melting) {
        /* Need T on element grid to do melting */
				p_to_centres(E,E->T,E->TE,E->mesh.levmax);

        /* Loop over elements to get melting by element */
				for(i=1;i<=E->mesh.nel;i++) {
					k = ((i-1) % E->mesh.elz) + 1;
					E->FmE[i] = melting(E,0.0,i,k);
				}

        /* Loop over nodes to reduce temperature */
				for(i=1;i<=E->mesh.nno;i++) {
					k = ((i-1) % E->mesh.noz) + 1;
					E->Fm[i] = melting_by_node(E,0.0,i,k);
				}

			  C = bulk_comp(E,4);

				/* Interpolate C,F from nodes to elements */
				/* Don't actually do this!
        p_to_centres(E,E->C,E->CE,E->mesh.levmax);
				p_to_centres(E,E->Fm,E->FmE,E->mesh.levmax);
 			  C = bulk_comp(E,41);
        */

				/* Update melt mass and bulk composition*/
				E->Total.bulk_comp_prev = E->Total.bulk_comp;
				E->Total.bulk_comp = 0.0;
			
				for(j=1;j<=E->mesh.esf;j++){ 
					E->slice.new_melt[j] = 0.0;

					for(k=1;k<=E->mesh.elz;k++) {
						e = (j-1)*E->mesh.elz + k;
						E->slice.melt[j] += E->FmE[e]*E->eco[e].area/E->Total.vol;
						E->slice.new_melt[j] += E->FmE[e]*E->eco[e].area;
						E->Total.melt_prod += E->FmE[e]*E->eco[e].area/E->Total.vol;
						E->Total.bulk_comp += E->CE[e] * E->eco[e].area;
/*						if (E->advection.timesteps == 470)
							fprintf(E->fpdebug,"a %d %g %g\n",e,E->CE[e],E->Total.bulk_comp);*/
					}

					E->slice.new_melt[j] /= E->Total.vol;
				}

        E->Total.bulk_comp /= E->Total.vol;

			}

			/* Impact Heating where appropriate */

			if (E->impacts.number >= 1)  
				for(i=0;i<E->impacts.number;i++) {
					if ((E->monitor.elapsed_time >= E->impacts.t[i]) && (E->impacts.flag[i] == 0)) {
						fprintf(E->fp,"Time for impact %d: %e\n",i,E->monitor.elapsed_time);  
						fprintf(stderr,"Time for impact %d: %e\n",i,E->monitor.elapsed_time);  
            fprintf(E->fp,"read heat from file? %d\n",E->impacts.heat_from_file);
            fprintf(stderr,"read heat from file? %d\n",E->impacts.heat_from_file);

            return_horiz_ave(E,E->T,E->Have.Tprev);
						impact_heating(E,i);

						E->impacts.flag[i]++; // Only impose heating once for each impact 
					}
				}

		}

    thermal_buoyancy(E);

    if( E->advection.timesteps < E->advection.max_timesteps)
        E->control.keep_going = 1;
    else
        E->control.keep_going = 0;

    on_off=(on_off==0)?1:0;

return;
}




void PG_timestep(E)
     struct All_variables *E;
{    
    void timestep();
    void predictor();
    void corrector();
    void pg_solver();
    void remove_horiz_ave();
    void return_horiz_ave();
    void std_timestep();
    void temperatures_conform_bcs();
		void impact_heating();
    void thermal_buoyancy();
    
    float Tmax(),T_interior1;

    int i,j,psc_pass,count,steps,iredo;
    int keep_going;

    float *DTdot, *T1, *Tdot1;

    static int loops_since_new_eta = 0;
    static int been_here = 0;
   
    DTdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    T1= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    Tdot1= (float *)malloc((E->mesh.nno+1)*sizeof(float));

    if (been_here++ ==0)    {
	  E->advection.timesteps=0;
      }

    E->advection.timesteps++;
 
    std_timestep(E);

    /* update temperature    */

       predictor(E,E->T,E->Tdot,1);
       for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++)   {
         pg_solver(E,E->T,E->Tdot,DTdot,E->V,1.0,E->TB,E->node);
         corrector(E,E->T,E->Tdot,DTdot,1);
         }	     

   if (E->control.composition)         {
	     predictor(E,E->C,E->Cdot,0);
	     for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++) {
	 	   pg_solver(E,E->C,E->Cdot,DTdot,E->V,E->control.comp_diff,E->CB,E->node);
		   corrector(E,E->C,E->Cdot,DTdot,0);
	       }	     
         }



    E->advection.total_timesteps++;

    temperatures_conform_bcs(E,E->T); 

		/* Apply impact heating if appropriate */
  	if (E->impacts.number >= 1) 
    	for(i=0;i<E->impacts.number;i++) {
      	if ((E->monitor.elapsed_time >= E->impacts.t[i]) && (E->impacts.flag[i] == 0)) {
          fprintf(E->fp,"Time for impact %d: %e\n",i,E->monitor.elapsed_time);  
          fprintf(stderr,"Time for impact %d: %e\n",i,E->monitor.elapsed_time); 
          fprintf(E->fp,"read heat from file? %d\n",E->impacts.heat_from_file);
          fprintf(stderr,"read heat from file? %d\n",E->impacts.heat_from_file);

          return_horiz_ave(E,E->T,E->Have.Tprev);
	        impact_heating(E,i);
  	      E->impacts.flag[i]++; // Only impose heating once for each impact 

    	  }
			}

    thermal_buoyancy(E);
 
    if( E->advection.timesteps < E->advection.max_timesteps)
        E->control.keep_going = 1;
    else
 	    E->control.keep_going = 0;

    free((void *) DTdot );  
    free((void *) T1 );   
    free((void *) Tdot1 );
    
return;  
}


/* ==============================
   predictor and corrector steps.
   ============================== */

void predictor(E,field,fielddot,ic)
     struct All_variables *E;
     float *field,*fielddot;
     int ic;

{ 
    int node;
    float multiplier;

   multiplier = (1.0-E->advection.gamma) * E->advection.timestep;

  if (ic==1)
    for(node=1;node<=E->mesh.nno;node++)  {
	if(!(E->node[node] & (OFFSIDE | TBX | TBZ | TBY))) 
	    field[node] += multiplier * fielddot[node] ;
	fielddot[node] = 0.0;
       }
  else
    for(node=1;node<=E->mesh.nno;node++)  {
	if(!(E->node[node] & OFFSIDE )) 
	    field[node] += multiplier * fielddot[node] ;
	fielddot[node] = 0.0;
       }

   return; }

void corrector(E,field,fielddot,Dfielddot,ic)
     struct All_variables *E;
     float *field,*fielddot,*Dfielddot;
     int ic;
    
{  int node;
   float multiplier;

   multiplier = E->advection.gamma * E->advection.timestep;

 if (ic==1)
   for(node=1;node<=E->mesh.nno;node++) {
       if(!(E->node[node] & (OFFSIDE | TBX | TBZ | TBY)))
	   field[node] += multiplier * Dfielddot[node];
       fielddot[node] +=  Dfielddot[node]; 
       }
 else
   for(node=1;node<=E->mesh.nno;node++) {
       if(!(E->node[node] & OFFSIDE ))
	   field[node] += multiplier * Dfielddot[node];
       fielddot[node] +=  Dfielddot[node]; 
       }
   
   return;  
 }

/* ===================================================
   The solution step -- determine residual vector from
   advective-diffusive terms and solve for delta Tdot
   Two versions are available -- one for Cray-style 
   vector optimizations etc and one optimized for 
   workstations.
   =================================================== */

//void pg_solver(E,T,Tdot,DTdot,V,diff,TBC,FLAGS)
void pg_solver(E,T,Tdot,DTdot,V,diff,TBC,FLAGS1) //
     struct All_variables *E;
     float *T,*Tdot,*DTdot;
     float **V;
     float diff;
     float **TBC;
     //unsigned int *FLAGS;
     unsigned int *FLAGS1; //
{
    void get_global_shape_fn();
    void pg_shape_fn();
    void element_residual();
    int el,e,a,i,a1;
    double xk[3][5],Eres[9];  /* correction to the (scalar) Tdot field */

    struct Shape_function PG;
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
 
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int ends=enodes[dims];

    
    for(i=1;i<=E->mesh.nno;i++)
 	  DTdot[i] = 0.0;

    for(el=1;el<=E->mesh.nel;el++)    {

	  i = (el-1)%E->mesh.elz+1;
	  /*diff =(E->diffusivity[i]+E->diffusivity[i+1])*0.5;*/

	  get_global_shape_fn(E,el,&GN,&GNx,&dOmega,xk,0,E->mesh.levmax);
	  pg_shape_fn(E,el,&PG,&GNx,V,diff);
	  element_residual(E,el,PG,GNx,dOmega,V,T,Tdot,Eres,diff);

      for(a=1;a<=ends;a++) {
	    a1 = E->ien[el].node[a];
	    DTdot[a1] += Eres[a]; 
        }

      } /* next element */

    for(i=1;i<=E->mesh.nno;i++) {
 	  if(E->node[i] & OFFSIDE) continue;
	  DTdot[i] *= E->Mass[i];         /* lumped mass matrix */
      }
     
    return;    
}



/* ===================================================
   Petrov-Galerkin shape functions for a given element
   =================================================== */

void pg_shape_fn(E,el,PG,GNx,V,diffusion)
     struct All_variables *E;
     int el;
     struct Shape_function *PG;
     struct Shape_function_dx *GNx;
     float **V;
     float diffusion;

{ 
    int i,j,node;
    int *ienmatrix;

    double uc1,uc2,uc3;
    double u1,u2,u3,VV[4][9];
    double uxse,ueta,ufai,xse,eta,fai,dx1,dx2,dx3,adiff;

    double prod1,unorm,twodiff;
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int lev=E->mesh.levmax;
    const int nno=E->mesh.nno;
    const int ends=enodes[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
  
    ienmatrix=E->ien[el].node;

    twodiff = 2.0*diffusion;
 
    uc1 =  uc2 = uc3 = 0.0;

    for(i=1;i<=ends;i++)   {
        node = ienmatrix[i];
        VV[1][i] = V[1][node];
        VV[2][i] = V[2][node];
      }
         
    for(i=1;i<=ENODES2D;i++) {
      uc1 +=  E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
      uc2 +=  E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
      }
    dx1 = 0.5*(E->X[1][ienmatrix[3]]+E->X[1][ienmatrix[4]]
              -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[2]]);
    dx2 = 0.5*(E->X[2][ienmatrix[3]]+E->X[2][ienmatrix[4]]
              -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[2]]);
    uxse = fabs(uc1*dx1*E->eco[el].centre[2]+uc2*dx2); 


    dx1 = 0.5*(E->X[1][ienmatrix[2]]+E->X[1][ienmatrix[3]]
              -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[4]]);
    dx2 = 0.5*(E->X[2][ienmatrix[2]]+E->X[2][ienmatrix[3]]
              -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[4]]);
    ueta = fabs(uc1*dx1*E->eco[el].centre[2]+uc2*dx2); 

    xse = (uxse>twodiff)? (1.0-twodiff/uxse):0.0;
    eta = (ueta>twodiff)? (1.0-twodiff/ueta):0.0;

    unorm = uc1*uc1 + uc2*uc2;

    adiff = (unorm>0.000001)?( (uxse*xse+ueta*eta)/(2.0*unorm) ):0.0;

    for(i=1;i<=VPOINTS2D;i++) {
          u1 = u2 = 0.0;
          for(j=1;j<=ENODES2D;j++)  /* this line heavily used */ {
            u1 += VV[1][j] * E->N.vpt[GNVINDEX(j,i)]; 
            u2 += VV[2][j] * E->N.vpt[GNVINDEX(j,i)];
	    }
	    
	  for(j=1;j<=ENODES2D;j++) {
             prod1 = (u1 * GNx->vpt[GNVXINDEX(0,j,i)] +
                      u2 * GNx->vpt[GNVXINDEX(1,j,i)]);
             PG->vpt[GNVINDEX(j,i)] = E->N.vpt[GNVINDEX(j,i)] + adiff * prod1;
	     }
	  }

   return;
 }



/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term.
   =========================================  */

void element_residual(E,el,PG,GNx,dOmega,V,field,fielddot,Eres,diff)
     struct All_variables *E;
     int el;
     struct Shape_function PG;
     struct Shape_function_dA dOmega;
     struct Shape_function_dx GNx;
     float **V;
     float *field,*fielddot;
     double Eres[9];
     float diff;

{
    int i,j,a,k,node,nodes[4],d,aid,back_front,onedfns;
    double Q;
    double dT[9],VV[4][9];
    double tx1[9],tx2[9],tx3[9];
    double v1[9],v2[9],v3[9];
    double adv_dT,t2[4];
    double T,DT;
    static int been_here=0;
 
    register double prod,sfn;
    struct Shape_function1 GM;
    struct Shape_function1_dA dGamma;
    double temp;

    void get_global_1d_shape_fn();
 
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    const int vpts=vpoints[dims];
    const int diffusion = (diff != 0.0); 

    for(i=1;i<=vpts;i++)	{ 
      dT[i]=0.0;
      v1[i] = tx1[i]=  0.0;
      v2[i] = tx2[i]=  0.0;
      }

	for(i=1;i<=ends;i++)   {
      node = E->ien[el].node[i];
        VV[1][i] = V[1][node];
        VV[2][i] = V[2][node];
      }
         
  
    for(j=1;j<=ends;j++)       {
      node = E->ien[el].node[j];
	  T = field[node];  
	  if(E->node[node] & (TBX | TBY | TBZ))
	    DT=0.0;
	  else
	    DT = fielddot[node];
	
	    for(i=1;i<=vpts;i++)  {
	 	  dT[i] += DT * E->N.vpt[GNVINDEX(j,i)];
		  tx1[i] +=  GNx.vpt[GNVXINDEX(0,j,i)] * T; 
		  tx2[i] +=  GNx.vpt[GNVXINDEX(1,j,i)] * T;   
		  sfn = E->N.vpt[GNVINDEX(j,i)];
		  v1[i] += VV[1][j] * sfn;
		  v2[i] += VV[2][j] * sfn;
	    }
      }


   Q=0;

/*
   if (diff>0.9)   {
     for(j=1;j<=ends;j++) 
       Q += E->C[E->ien[i].node[j]]; 
     Q = Q/ends;
     Q = Q*E->control.Q0;
     }
*/

/* Only put rad_heat in silicate layer, not ice */
   if (E->eco[el].centre[2] <= E->viscosity.zlm)
     Q = E->rad_heat.total;
   else
     Q = 0.0;

/* Add rad_heat everywhere */
// Q = E->rad_heat.total;

/*
   if (E->control.composition)
	    Q = (1-E->CE[el])*E->rad_heat.total+E->CE[el]*E->control.Q0ER;
*/
   if(E->control.tidal_heating==0) {
     Q = (Q + E->heating_visc[el] - E->heating_adi[el])/E->heating_latent[el];
     }
   else {

     Q = (Q + E->heating_visc[el] - E->heating_adi[el] + E->heating_tidal[el]*E->tidal_visc[el])
	/ E->heating_latent[el];
     }

   if(E->control.shear_heating)
		Q += E->heating_shear[el] / E->heating_latent[el];

   if(E->control.despin && !(E->control.despun))
		Q += E->heating_despin[el]*E->tidal_visc[el] / E->heating_latent[el];

    /* construct residual from this information */


    if(diffusion){
       for(j=1;j<=ends;j++) {
          Eres[j]=0.0;
          for(i=1;i<=vpts;i++) 
            Eres[j] -= PG.vpt[GNVINDEX(j,i)] * dOmega.vpt[i] 
	             * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i]) +
                     diff/E->heating_latent[el]*dOmega.vpt[i] 
                     * (GNx.vpt[GNVXINDEX(0,j,i)] * tx1[i] +
		        GNx.vpt[GNVXINDEX(1,j,i)] * tx2[i] ); 

	    }
      }

    else { /* no diffusion term */
	    for(j=1;j<=ends;j++) {
    	  Eres[j]=0.0;
		  for(i=1;i<=vpts;i++) 
		    Eres[j] -= PG.vpt[GNVINDEX(j,i)] * dOmega.vpt[i] * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i]);
		  }
      }	

	/* See brooks etc: the diffusive term is excused upwinding for 
	   rectangular elements  */
  
    /* include BC's for fluxes at (nominally horizontal) edges (X-Y plane) */

/*    if(FLAGS!=NULL) {
	onedfns=0;
	for(a=1;a<=ends;a++)
	    if (FLAGS[E->ien[el].node[a]] & FBZ) {
		if (!onedfns++) get_global_1d_shape_fn(E,el,&GM,&dGamma);
 
		nodes[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
		nodes[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
		nodes[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][2];
		nodes[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][2];
	  
		for(aid=0,j=1;j<=onedvpoints[E->mesh.nsd];j++)
		    if (a==nodes[j])
			aid = j;
		if(aid==0)  
		    printf("%d: mixed up in pg-flux int: looking for %d\n",el,a);

		if (loc[a].plus[1] != 0)
		    back_front = 0;
		else back_front = dims;

		for(j=1;j<=onedvpoints[dims];j++)
		    for(k=1;k<=onedvpoints[dims];k++)
			Eres[a] += dGamma.vpt[GMVGAMMA(1+back_front,j)] *
			    E->M.vpt[GMVINDEX(aid,j)] * g_1d[j].weight[dims-1] *
			    BC[2][E->ien[el].node[a]] * E->M.vpt[GMVINDEX(k,j)];
	    }
    } 
    
 */
    return; 
}




/* =====================================================
   Obtain largest possible timestep (no melt considered)
   =====================================================  */


void std_timestep(E)
     struct All_variables *E;

{ 
    static int been_here = 0;
    static float diff_timestep,root3,root2;
    int i,d,n,nel,el,node;

    float adv_timestep;
    float ts,uc1,uc2,uc3,uc,size,step,VV[4][9];
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    
	nel=E->mesh.nel;

    if(E->advection.fixed_timestep != 0.0) {
      E->advection.timestep = E->advection.fixed_timestep;
      return;
    }

    if (been_here == 0)  {
	  diff_timestep = 1.0e8; 
	  for(el=1;el<=nel;el++)  { 
	      ts = E->eco[el].size[1]*E->eco[el].size[1]*E->eco[el].centre[2]*E->eco[el].centre[2];
	      diff_timestep = min(diff_timestep,ts);
	      ts = E->eco[el].size[2]*E->eco[el].size[2];
	      diff_timestep = min(diff_timestep,ts);
	    }
      diff_timestep = 0.5 * diff_timestep;
      }


    
  adv_timestep = 1.0e8;
  for(el=1;el<=nel;el++) {

	for(i=1;i<=ends;i++)   {
      node = E->ien[el].node[i];
        VV[1][i] = E->V[1][node];
        VV[2][i] = E->V[2][node];
        if(dims==3) VV[3][i] = E->V[3][node];
      }
         
    uc=uc1=uc2=uc3=0.0;	  
    if(3==dims) {
      for(i=1;i<=ENODES3D;i++) {
        uc1 += E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
        uc2 += E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
        uc3 += E->N.ppt[GNPINDEX(i,1)]*VV[3][i];
        }
      uc = fabs(uc1)/E->eco[el].size[1] + fabs(uc2)/E->eco[el].size[2] + fabs(uc3)/E->eco[el].size[3];

	  step = (0.5/uc);
	  adv_timestep = min(adv_timestep,step);
      }
    else {
	  for(i=1;i<=ENODES2D;i++) {
        uc1 += E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
        uc2 += E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
        }
      uc = fabs(uc1)/(E->eco[el].size[1]*E->eco[el].centre[2]) + fabs(uc2)/E->eco[el].size[2];
	  
	  step = (0.5/uc);
	  adv_timestep = min(adv_timestep,step);
      }
    }

    adv_timestep = E->advection.dt_reduced * adv_timestep;         

    adv_timestep =  1.0e-32+min(E->advection.fine_tune_dt*adv_timestep,diff_timestep);
    
    E->advection.timestep = adv_timestep;

    return; 
  }


  void process_heating(E)
  struct All_variables *E;
  { 

  int e,i,j,ee;
  static int been=0;
  static float para1;
  double relem,thelem,x;
  double temp1,temp2,temp3;
  double *viscr;
  FILE *fp;
  char filename[250];
  void return_horiz_ave();


    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->mesh.nno;
    const int vpts = vpoints[dims];
    const double two=2.0;

  void strain_rate_2_inv();
	void despin_heat();

  viscr = (double *) malloc((E->mesh.nel+1)*sizeof(double));

  if (been==0)  {
    para1 = E->sphere.ro_dim*E->sphere.ro_dim/(E->data.density*E->data.Cp*E->data.ref_temperature*E->data.therm_diff);
//    E->data.disptn_number = E->data.therm_exp*E->data.grav_acc*E->sphere.ro_dim/E->data.Cp;
    been++;
    }

  return_horiz_ave(E,E->T,E->Have.T);

  if (E->rad_heat.num==0)  {
     E->rad_heat.total = E->control.Q0; 
  }
  else if (E->rad_heat.num==1)  {
     E->rad_heat.total = para1*E->control.Q0; 
  }
  else if (E->rad_heat.num==2) {
     temp1 = 4.6e9-E->monitor.time_scale*E->monitor.elapsed_time;
     temp2 = E->rad_heat.percent[0]*E->rad_heat.concen[0]*E->rad_heat.heat_g[0]
	         *exp(temp1*log(two)/E->rad_heat.decay_t[0])
           + E->rad_heat.percent[1]*E->rad_heat.concen[1]*E->rad_heat.heat_g[1]
	         *exp(temp1*log(two)/E->rad_heat.decay_t[1])
           + E->rad_heat.percent[2]*E->rad_heat.concen[2]*E->rad_heat.heat_g[2]
	         *exp(temp1*log(two)/E->rad_heat.decay_t[2])
           + E->rad_heat.percent[3]*E->rad_heat.concen[3]*E->rad_heat.heat_g[3]
	         *exp(temp1*log(two)/E->rad_heat.decay_t[3]);
     E->rad_heat.total = para1*E->data.density*temp2; 
  }

  temp1 = E->data.disptn_number/E->control.Ra_temp;

  for (e=1;e<=E->mesh.nel;e++)  {
    E->heating_latent[e] = 1.0;
    }

 if (E->control.visc_heating)  {
   strain_rate_2_inv(E,E->heating_visc,0);
 
   for (e=1;e<=E->mesh.nel;e++)  {
     temp2 = 0.0;
     for (i=1;i<=vpts;i++)
       temp2 += E->EVi[(e-1)*vpts+i];
     temp2 = temp2/vpts;

     E->heating_visc[e] = temp1*temp2*E->heating_visc[e];
     }
   }

 if (E->control.adi_heating)  {
   for (e=1;e<=E->mesh.nel;e++)  {
    ee = (e-1)%E->mesh.elz+1;

    temp2 = 0.0;
    for (i=1;i<=ends;i++)    {
      j = E->ien[e].node[i];
      temp2 = temp2 + E->V[2][j]*(E->T[j]+E->data.surf_temp )*E->data.disptn_number;  
      }
    temp2 = temp2/ends;
    E->heating_adi[e] = temp2*(E->expansivity[ee]+E->expansivity[ee+1])*0.5;

    }
  }

/* 
 * The tidal heating model used to generate the initial heating read in by
 * Citcom assumes constant material properties.  The viscosity varies
 * laterally and may enhance or reduce the heating in a particular element.
 * The following statement accounts for that (Ojakangas and Stevenson, 1989).
 *
 * This is really only designed to work off of ref_visc as the viscosity used
 * for all layers in the tidal code. If you have many layers with different 
 * viscosities, you'll need to modify this to use a reference viscosity for 
 * each layer.
 */

  if (E->control.tidal_heating || E->control.despin) {

    for (e=1;e<=E->mesh.nel;e++)  {
      ee = (e-1)%E->mesh.elz+1;
//      relem = 0.5*(E->X[2][ee] + E->X[2][ee+1]);
      relem = E->eco[e].centre[2];

      /* Find average viscosity at radial position of e *
       * Important quantity is lateral variation of viscosity
       * Use log average due to large variations
       */
      viscr[ee] = exp( 0.5*( log(E->Have.Vi[ee]) + log(E->Have.Vi[ee+1]) ) );
      
      E->tidal_visc[e] = 1.0;
      temp3 = E->data.frequency * E->data.ref_viscosity / E->data.rigidity;
      temp3 = temp3*viscr[ee];

      temp2 = 0.0;
      for (i=1;i<=vpts;i++)
        temp2 += E->EVi[(e-1)*vpts+i];
      temp2 = temp2/vpts;

      temp2 = temp2/viscr[ee];  // Care about lateral variation

      if (relem > E->viscosity.zlith) /* Modify for lithosphere */
	temp1 = E->viscosity.max_value;
      else 
	temp1 = 1.0;

        E->tidal_visc[e] = 1.0; /* Use this only if convection is weak */
      E->tidal_visc[e] = temp2 * (1 + temp3*temp3) / (1 + temp2*temp2*temp3*temp3);

/*      E->tidal_visc[e] = (temp2/temp1) * (1 + temp1*temp1*temp3*temp3) 
			/ (1 + temp2*temp2*temp3*temp3);
*/
//    if (e<=65)
//      fprintf(stderr,"tidal visc: %d %g %g %g %g\n",e,E->data.frequency,E->data.ref_viscosity,E->data.rigidity,E->tidal_visc[e]);
//      fprintf(stderr,"tidal visc: %d %d %g %g %g\n",e,ee,temp2,temp3,E->tidal_visc[e]);
      
    }
  }

  if(E->control.despin && !(E->control.despun))
		despin_heat(E);


  if (E->control.Ra_670!=0.0)   {
    temp1 = 2.0*E->control.clapeyron670*E->control.Ra_670/(E->control.Ra_temp/E->control.width670);
    for (e=1;e<=E->mesh.nel;e++)  {
      temp2 = 0;
      temp3 = 0;
      for (i=1;i<=ends;i++)    {
        j = E->ien[e].node[i];
        temp2 = temp2 + temp1*(1.0-E->Fas670[j])*E->Fas670[j]
                   *E->V[2][j]*(E->T[j]+E->data.surf_temp)*E->data.disptn_number;
        temp3 = temp3 + temp1*E->control.clapeyron670
                        *(1.0-E->Fas670[j])*E->Fas670[j]
                        *(E->T[j]+E->data.surf_temp)*E->data.disptn_number;
        }
      temp2 = temp2/ends;
      temp3 = temp3/ends;
      E->heating_adi[e] += temp2;
      E->heating_latent[e] += temp3;
      }
    }

  if (E->control.Ra_410!=0.0)   {
    temp1 = 2.0*E->control.clapeyron410*E->control.Ra_410/(E->control.Ra_temp/E->control.width410);
    for (e=1;e<=E->mesh.nel;e++)  {
      temp2 = 0;
      temp3 = 0;
      for (i=1;i<=ends;i++)    {
        j = E->ien[e].node[i];
        temp2 = temp2 + temp1*(1.0-E->Fas410[j])*E->Fas410[j]
                   *E->V[2][j]*(E->T[j]+E->data.surf_temp)*E->data.disptn_number;
        temp3 = temp3 + temp1*E->control.clapeyron410
                       *(1.0-E->Fas410[j])*E->Fas410[j]
                       *(E->T[j]+E->data.surf_temp)*E->data.disptn_number;
        }
      temp2 = temp2/ends;
      temp3 = temp3/ends;
      E->heating_adi[e] += temp2;
      E->heating_latent[e] += temp3;
      }
    }

  fprintf(E->fp,"QQ %g %g %g\n",E->rad_heat.total,E->data.crust_rad*E->rad_heat.total,E->data.mantle_rad*E->rad_heat.total);

if((E->monitor.solution_cycles%E->control.record_every)==0)  {
  sprintf(filename,"%s/heating.%d",E->control.data_file,E->monitor.solution_cycles);  
  fp=fopen(filename,"w");
  fprintf(fp,"QQ %g %g %g\n",E->control.Ra_temp,E->data.disptn_number,E->rad_heat.total);
  for (e=1;e<=E->mesh.nel;e++){
    ee = (e-1)%E->mesh.elz+1;
    fprintf(fp,"%d %g %g %g %g %g %g %d %g\n",e,E->EVi[(e-1)*vpts+1],viscr[ee],E->heating_visc[e],E->heating_adi[e],E->heating_tidal[e]*E->tidal_visc[e],E->tidal_visc[e],(E->eco[e].centre[2] <= E->viscosity.zlith),E->heating_tidal[e]);
/*    fprintf(fp,"%d %g %g %g %g %g %g\n",e,E->EVi[(e-1)*vpts+1],
			E->heating_visc[e],E->heating_adi[e],
			(E->heating_tidal[e]*E->tidal_visc[e]),
			(E->heating_despin[e]*E->tidal_visc[e]),E->tidal_visc[e]);*/
      }
 fclose(fp);
   }
/*
*/

  free((void *) viscr );

  return; 
  }

/****************************************************************
 * despin_heat																									*
 * 																															*
 * Computes the stresses and heating due to despinning of a 		*
 * body based on the orbital evolution and mechanical 					*
 * properties.																									*
 *																															*
 * This is for only the current timestep, not an average				*
 ****************************************************************/

void despin_heat(E)
  struct All_variables *E;
  {
		
  FILE *fp1; 
  char output_file[255];
	int e,i,j,k,n;
	double m1,m2, rot_prev,rot;
	long double semimajor,semimajor_prev,mean_motion;
	long double Ls1,Ls2,Lo1,Lo2,DLs;
	long double Es1,Es2,Eo1,Eo2;
	double Ediss,dt;
	float poisson;
	double mu,omega,eta,drotdt;
	complex double c_rigidity, c_k2;
	complex double f1,f2;
	complex double *log_c_Ri;
	complex double ave_c_Ri;
	double volume, e_volume;

  int nn[7];
  float coeff_beta[7];
	float density_ratio;
	complex double kappa_grav;
	complex double A1,A2,B1,B2,D1,D2,E1,E2;
	complex double *M_mat[5], *M_inv[5];
  complex double K_vec[5], weight[7];
  complex double u_r, u_th, du_r_dr, du_th_dth, du_th_dr, du_r_dth;
  complex double e_rr, e_tt, e_ff, e_rt; 
  complex double s_rr, s_tt, s_ff, s_rt; 
  double relem,thelem,x;
  double *work, global_diss;

	void gaussj();

  log_c_Ri = (complex double *) malloc((E->mesh.nno+1)*sizeof(complex double));
  work = (double *) malloc((E->mesh.nel+1)*sizeof(double));
  for(i=1;i<=4;i++) {
		M_mat[i] = (complex double *)  malloc(5.0*sizeof(complex double));
		M_inv[i] = (complex double *)  malloc(5.0*sizeof(complex double));
	}

  global_diss = 0.0;
	poisson = 0.5;	/* Assume incompressibility */
	if (poisson = 0.5)
		poisson -= 1.0e-6;	/* Avoid div. by zero */
	density_ratio = E->data.density / E->data.density_core;

	/* Assemble complex rigidity array based on viscosity */
		
	omega = E->data.rotation;
	mu = E->data.rigidity;
  for (n=1;n<=E->mesh.nno;n++)  {	/* By Node */
		eta = E->Vi[n]*E->data.ref_viscosity;
 		c_rigidity = (omega*eta + I*mu) ;
 		c_rigidity *= (mu*omega*eta)/(mu*mu + omega*omega*eta*eta);
		E->c_Ri[n] = c_rigidity;
	}
  for (e=1;e<=E->mesh.nel;e++)  {	/* By Element */
		eta = E->EVi[e]*E->data.ref_viscosity;
  	c_rigidity = (omega*eta + I*mu) ;
  	c_rigidity *= (mu*omega*eta)/(mu*mu + omega*omega*eta*eta);
		E->c_ERi[e] = c_rigidity;
	}

	/* Get some kind of average */
	
  for (n=1;n<=E->mesh.nno;n++){
    log_c_Ri[n] = clog(E->c_Ri[n]);
  }

	/* Use layers when we get a prop. mat. soln in place for k2.  For now
   * just take radial ave. of rig. as well, weighted by Mass Matrix */

	/*
  return_horiz_ave(E,log_r_Ri,E->Have.r_Ri);
  return_horiz_ave(E,log_i_Ri,E->Have.i_Ri);
  for (k=1;k<=E->mesh.noz;k++) {
		E->Have.r_Ri[k] = exp(E->Have.r_Ri[k]);
		E->Have.i_Ri[k] = exp(E->Have.i_Ri[k]);
	}
	 */

	ave_c_Ri = 0.0;
	volume = 0.0;
/*	(void)fprintf(E->fp,"n %d, volume %g\n",n,volume);*/

	for (n=1;n<=E->mesh.nno;n++) {
		ave_c_Ri += log_c_Ri[n] / E->Mass[n]; /* Really Ri * dVol. */
		volume += 1.0/E->Mass[n];
	/*	(void)fprintf(E->fp,"n %d, Mass %g, volume %g, log_c_Ri %g + %g i, ave_c_Ri %g + %gi\n",
				n,E->Mass[n],volume,creal(log_c_Ri[n]),cimag(log_c_Ri[n]),
				creal(ave_c_Ri),cimag(ave_c_Ri));*/
	}
	ave_c_Ri = ave_c_Ri/volume;
	c_rigidity = cexp(ave_c_Ri);

  c_k2 = 1.5 / (1.0 + (19.0) * c_rigidity 
											/ (2*E->data.grav_acc*E->data.density*E->sphere.ro_dim));
  drotdt = (3.0 * -cimag(c_k2) * E->data.grav_const * E->data.mass_primary 
							* E->data.mass_primary * pow((E->sphere.ro_dim),5)) 
							/ (E->data.moi * pow(E->data.semimajor_axis,6));

	(void)fprintf(E->fp,"ave_c_Ri %g + %g i\n rig %g + %g i, k2 %g + %g i\n drotdt %g\n",	
			creal(ave_c_Ri),cimag(ave_c_Ri),
			creal(c_rigidity),cimag(c_rigidity),creal(c_k2),cimag(c_k2),
			drotdt);

	/* New rotation rate */
	rot_prev = E->data.rotation;
	dt = E->advection.timestep 
				* (E->sphere.ro_dim * E->sphere.ro_dim / E->data.therm_diff);
	rot = rot_prev - E->data.despin_rate * dt;

  if(rot <= E->data.rot_final) {
		rot = E->data.rot_final;
		E->control.despun = 1;
		E->control.despin_timescale = E->monitor.elapsed_time
									* (E->sphere.ro_dim * E->sphere.ro_dim / E->data.therm_diff);
	}

	E->data.rotation = rot;

	(void)fprintf(E->fp,"rot_prev %.15e -> rot %.15e, drot = %.9e, despin =%e, dt=%e -> %.10e\n",
			rot_prev,rot,rot-rot_prev,E->data.despin_rate,E->advection.timestep,dt);

	/* 
	 * Compute global dissipation rate at this timestep
	 * We know the moment of inertia, initial and final rotation 	
	 * rates, mass, and final orbital semimajor axis. 
	 */

	 /* Solve for new semimajor axis using conservation of angular 	momentum */

	semimajor_prev = E->data.semimajor_axis;
  
	Ls1 = E->data.moi*rot_prev;
  Ls2 = E->data.moi*rot;
	Lo1 = E->data.mass*E->data.mean_motion*semimajor_prev*semimajor_prev;
	DLs = Ls2-Ls1;
	
  Lo2 = Lo1 - DLs;
/*	fprintf(stderr,"Ls1 %.16Lg Ls2 %.16Lg \nLo1 %.16Lg Lo2 %.16Lg\n",Ls1,Ls2,Lo1,Lo2);
	fprintf(stderr,"DLs %Lg DLo %Lg DL %Lg\n",DLs,Lo2-Lo1,DLs+Lo2-Lo1);
	fprintf(stderr,"Mass %g rot_prev %g a_prev %Lg\n",E->data.mass,rot_prev,semimajor_prev);
*/
	if(E->advection.timestep == 0.0)
		semimajor = semimajor_prev;
	else
	  semimajor = powl((Lo2/E->data.mass),2) 
		 											/ (E->data.grav_const * E->data.mass_primary);

	E->data.semimajor_axis = semimajor;
  E->data.mean_motion = sqrt((E->data.grav_const * E->data.mass_primary) / 
								  pow(semimajor,3));

	(void)fprintf(stderr,"a_prev %.12Lg -> a %.12Lg, da = %Lg\n",
			semimajor_prev, semimajor, semimajor-semimajor_prev);

  /* Next, solve for mean motion using Kepler's 3rd Law */
 	mean_motion = sqrt((E->data.grav_const * E->data.mass_primary) / 
											  pow(semimajor,3));

 	/* Find dissipated energy using conservation of energy */
	Es1 =  E->data.moi*rot_prev*rot_prev / 2.0;
	Es2 =  E->data.moi*rot*rot / 2.0;
	Eo1 = - E->data.grav_const*E->data.mass_primary*E->data.mass 
				/ (2.0*semimajor_prev);
	Eo2 = - E->data.grav_const*E->data.mass_primary*E->data.mass 
				/ (2.0*semimajor);
/*	fprintf(stderr,"Es1 %Lg Es2 %Lg \nEo1 %.16Lg Eo2 %.16Lg\n",Es1,Es2,Eo1,Eo2);
	fprintf(stderr,"Es1 - Es2 %Lg Eo1 - Eo2 %Lg\n",Es1-Es2,Eo1-Eo2);
*/
	Ediss = ( E->data.moi*(rot_prev*rot_prev - rot*rot) 
			  				- E->data.grav_const*E->data.mass_primary*E->data.mass 
		  							* (1.0/semimajor_prev - 1.0/semimajor) )/2.0;

	E->data.Ediss += Ediss;  /* Keep track of dissipated energy */	
  fprintf(stderr,"Ediss = %g, Total = %g, dt=%e\n",Ediss,E->data.Ediss,dt);
  fprintf(E->fp,"Ediss = %g, Total = %g, dt=%e\n",Ediss,E->data.Ediss,dt);

  /* Flattening */
  m1 = (rot_prev*rot_prev*E->sphere.ro_dim / E->data.grav_acc);
  m2 = (rot*rot*E->sphere.ro_dim / E->data.grav_acc);
  f1 = E->data.flattening;
  f2 = f1 - (1.25*(m1-m2)) / (1.0 + (19.0*E->data.rigidity) 
						/(2.0*E->data.grav_acc*E->data.density*E->sphere.ro_dim));
  fprintf(E->fp,"m = %g -> %g\n",m1,m2);

  /* Solution parameters */
	kappa_grav = (4.0*M_PI/3.0) * E->data.grav_const * E->sphere.ro_dim 
		* E->sphere.ro_dim * E->data.density * E->data.density / c_rigidity;

  fprintf(E->fp,"kappa_grav = %.16g + %.16g\n",
			creal(kappa_grav),cimag(kappa_grav));
  nn[1] = 1;
  nn[2] = 3;
  nn[3] = -2;
  nn[4] = -4;

	/* Assuming incompressibility */
  coeff_beta[1] = -3.0;
  coeff_beta[2] = (-7.0+4.0*poisson)/(2.0*poisson);
  coeff_beta[3] = -6.0*(1.0-2.0*poisson)/(5.0-4.0*poisson);
  coeff_beta[4] = 2.0;
  
	K_vec[1] = (kappa_grav / 12.0) * (m1 - m2);
	K_vec[2] = (kappa_grav / 12.0) * (m1 - m2) * (1.0 - density_ratio);
	K_vec[3] = 0.0;
	K_vec[4] = 0.0;
/*
	for(i=1;i<=4;i++) {
		fprintf(E->fp,"n %d beta %g K %g + %g i\n",nn[i],coeff_beta[i],
				creal(K_vec[i]),cimag(K_vec[i]));
	}
*/
	A1 = (0.1 - (1.0 + (density_ratio - 1.0) * pow(E->sphere.ri,3.0)) / 6.0) 
				* kappa_grav;
	B1 = 0.1 * (density_ratio - 1.0) * kappa_grav;
	D1 = 0.0;
	E1 = 0.0;
	A2 = 0.1 * (1.0 - density_ratio) * pow(E->sphere.ri,2.0) * kappa_grav;
	B2 = (density_ratio - 1.0) * (1.5 + density_ratio) * kappa_grav / 30.0;
	D2 = 0.0;
	E2 = 0.0;
/*
	fprintf(E->fp,"A1 %g + %g i A2 %g + %g i\n",
			creal(A1),cimag(A1),creal(A2),cimag(A2));
	fprintf(E->fp,"B1 %g + %g i B2 %g + %g i\n",
			creal(B1),cimag(B1),creal(B2),cimag(B2));
	fprintf(E->fp,"D1 %g + %g i D2 %g + %g i\n",
			creal(D1),cimag(D1),creal(D2),cimag(D2));
	fprintf(E->fp,"E1 %g + %g i E2 %g + %g i\n",
			creal(E1),cimag(E1),creal(E2),cimag(E2));
*/
	M_mat[1][1] = -2.0 + 6.0*(A1 + B1*pow(E->sphere.ri,5.0));
	M_mat[1][2] = 1.0 + 6.0*(A1 + B1*pow(E->sphere.ri,7.0));
	M_mat[1][3] = 6.0 + 6.0*(A1 + B1*pow(E->sphere.ri,2.0));
	M_mat[1][4] = 8.0 + 6.0*(A1 + B1);

	M_mat[2][1] = -2.0 + 6.0*(A2 + B2*pow(E->sphere.ri,2.0));
	M_mat[2][2] = pow(E->sphere.ri,2.0) + 6.0*(A2 + B2*pow(E->sphere.ri,4.0));
	M_mat[2][3] = 6.0*pow(E->sphere.ri,-3.0) + 6.0*(A2 + B2/E->sphere.ri);
  M_mat[2][4] = 8.0*pow(E->sphere.ri,-5.0) 
									+ 6.0*(A2 + B2*pow(E->sphere.ri,-3.0));

	M_mat[3][1] = 6.0 + 6.0*(D1 + E1*pow(E->sphere.ri,5.0));
	M_mat[3][2] = 16.0 + 6.0*(D1 + E1*pow(E->sphere.ri,7.0));
	M_mat[3][3] = 6.0 + 6.0*(D1 + E1*pow(E->sphere.ri,2.0));
	M_mat[3][4] = 16.0 + 6.0*(D1 + E1);
	
	M_mat[4][1] = 6.0 + 6.0*(D2 + E2*pow(E->sphere.ri,2.0));
	M_mat[4][2] = 16.0*pow(E->sphere.ri,2.0) 
									+ 6.0*(D2 + E2*pow(E->sphere.ri,4.0));
	M_mat[4][3] = 6.0*pow(E->sphere.ri,-3.0) + 6.0*(D2 + E2/E->sphere.ri);
  M_mat[4][4] = 16.0*pow(E->sphere.ri,-5.0) 
									+ 6.0*(D2 + E2*pow(E->sphere.ri,-3.0));
/*
	for(i=1;i<=4;i++) {
		for(j=1;j<=4;j++)
			fprintf(E->fp,"M %g + %g i\n",creal(M_mat[i][j]),cimag(M_mat[i][j]));
	}
*/
	/* For homog. planet */
	/*weight[2] = -((E->data.density*E->data.grav_acc)
				/(95.0 * c_rigidity * E->sphere.ro_dim)) * (f2 - 1.25*m2);
  weight[1] = -(8.0/3.0) * weight[2] * pow(E->sphere.ro_dim,2);

  fprintf(E->fp,"weight = %g + %g i, %g + %g i\n",
					creal(weight[1]),cimag(weight[1]),creal(weight[2]),cimag(weight[2]));
*/
	/* For layered planet, solve M_mat*weight_vec = K_vec */

	for(i=1;i<=4;i++) {
		for(j=1;j<=4;j++){ 
			M_inv[i][j] = M_mat[i][j];
		/*	M_inv[i][j] = 1.0;	// Test gaussj */
		}
		weight[i] = K_vec[i];
	}
/*
	M_inv[2][2] = 2.0;
	M_inv[2][4] = 2.0;
	M_inv[3][4] = 0.0;
	M_inv[4][2] = 4.0;
	M_inv[4][3] = 2.0;
	M_inv[4][4] = 3.0;

	weight[1] = 4.0;		
	weight[2] = 6.0;		
	weight[3] = 3.0;		
	weight[4] = 10.0;		
	*/
/*
	fprintf(E->fp,"M =\n");
  for(i=1;i<=4;i++) {
  	for(j=1;j<=4;j++) {
			fprintf(E->fp,"%g+%gi \t",creal(M_inv[i][j]),cimag(M_inv[i][j]));
		}
		fprintf(E->fp,"%g+%gi\n",creal(weight[i]),cimag(weight[i]));
	}
*/
	gaussj(E,M_inv,4,weight);
/*
	fprintf(E->fp,"inv(M) =\n");
  for(i=1;i<=4;i++) {
  	for(j=1;j<=4;j++) {
			fprintf(E->fp,"%g+%gi \t",creal(M_inv[i][j]),cimag(M_inv[i][j]));
		}
		fprintf(E->fp,"%g+%gi\n",creal(weight[i]),cimag(weight[i]));
	}
*/
	/* Normalize weighting by appropriate powers of radius, calculate new 
	 	 flattening */
	f2 = f1;
	for(i=1;i<=4;i++) {
		f2 += 6*weight[i];
		weight[i] *= pow(E->sphere.ro_dim,(1-nn[i]));	
	/*	fprintf(E->fp,"%g+%gi\n",creal(weight[i]),cimag(weight[i]));*/
		}
  fprintf(E->fp,"f = %g+%gi -> %g+%gi, df = %g+%gi\n",creal(f1),cimag(f1),
			creal(f2),cimag(f2),creal(f2-f1),cimag(f2-f1));
	
	/* quantities by element */
/*
  sprintf(output_file,"%s/despin.%d",E->control.data_file,
			E->monitor.solution_cycles);
  fp1=fopen(output_file,"w");
  fprintf(fp1,"Despin Quantities\n");
*/
	e_volume = 0.0;
  for (e=1;e<=E->mesh.nel;e++)  {
		work[e] = 0.0;
		/* E->heating_despin[e] = E->power_despin_ave; */ /* Average Value */
		/* Perturb based on stress pattern */
	
    relem = E->eco[e].centre[2]*E->sphere.ro_dim;
		thelem =  E->eco[e].centre[1];
		x = 1.0 + 3.0 * cos(2.0*thelem);
	
		/*if (e==1) 
    	fprintf(stderr,"relem = %g, thelem = %g, x = %g\n",relem,thelem,x);*/
		
		/* Displacement */
		u_r = 0.0;
		u_th = 0.0;
		du_r_dr = 0.0;
		du_th_dth = 0.0;
		du_th_dr = 0.0;
		du_r_dth = 0.0;

		for(i=1; i<=4; i++) {
  		u_r += weight[i] * pow(relem,nn[i]) * x;
		  u_th += weight[i] * coeff_beta[i] * pow(relem,nn[i]) * sin(2.0*thelem);
		  du_r_dr += weight[i] * nn[i] * pow(relem,(nn[i]-1)) * x;
  		du_th_dth += weight[i] * coeff_beta[i] * pow(relem,nn[i]) 
															* 2.0*cos(2.0*thelem);
		  du_th_dr += weight[i] * coeff_beta[i] * nn[i] * pow(relem,(nn[i]-1)) 
														* sin(2.0*thelem);
		  du_r_dth += weight[i] * pow(relem,nn[i]) * -6*sin(2.0*thelem);
		}

		/* Strain */
		e_rr = du_r_dr;
		e_tt = (du_th_dth + u_r) / relem;
		if((E->eco[e].centre[1] == 0.0) || (E->eco[e].centre[1] == M_PI))
  		e_ff = 0.0;
    else
		  e_ff = (u_th/tan(thelem) + u_r) / relem;
		e_rt = du_th_dr + (du_r_dth - u_th) / relem;

		/* Stress */
		s_rr = 0.0;
		s_tt = 0.0;
		s_ff = 0.0;
		s_rt = 0.0;

		for(i=1; i<=4; i++) {
			s_rr += pow(relem,(nn[i]-1)) * ( ((1.0 - poisson)*nn[i] + 2.0*poisson) 
																		+ poisson * coeff_beta[i]) * weight[i];
			s_tt += pow(relem,(nn[i]-1)) * ( ((1.0 + nn[i]*poisson) 
																					+ poisson*coeff_beta[i]) 
																			- (3.0*(1.0 + nn[i]*poisson) 
																					+ (2.0 - poisson) * coeff_beta[i]) 
																			* -cos(2*thelem) ) * weight[i];
			s_ff += pow(relem,(nn[i]-1)) * ( ((1.0 + nn[i]*poisson) 
																					+ (1.0 - poisson)*coeff_beta[i]) 
																			- (3.0*(1.0 + nn[i]*poisson) 
																					+ (1.0 + poisson) * coeff_beta[i]) 
																			* -cos(2*thelem) ) * weight[i];
			s_rt += pow(relem,(nn[i]-1)) * ((nn[i]-1)*coeff_beta[i]/2.0 - 3.0) 
								* weight[i];
		}

		s_rr *= 2*c_rigidity/(1-2*poisson) * x;
		s_tt *= 2*c_rigidity/(1-2*poisson);
		s_ff *= 2*c_rigidity/(1-2*poisson);
		s_rt *= 2*c_rigidity * sin(2*thelem);


		/* This is really only good for diagonal term */
		/*
		s_rr = 2.0*c_rigidity*e_rr;
		s_tt = 2.0*c_rigidity*e_tt;
		s_ff = 2.0*c_rigidity*e_ff;
		s_rt = 2.0*c_rigidity*e_rt;
		*/
/*
		if(e==1) {
		  		fprintf(stderr,"ur = %g + %gi, ut = %g + %gi\n",creal(u_r),cimag(u_r),
					creal(u_th),cimag(u_th));
				  fprintf(stderr,"durdr = %g + %gi, dutdt = %g + %gi, dutdr = %g + %gi, durdt = %g + %gi \n",
					creal(du_r_dr),cimag(du_r_dr), creal(du_th_dth),cimag(du_th_dth),
					creal(du_th_dr),cimag(du_th_dr), creal(du_r_dth),cimag(du_r_dth));
				  fprintf(stderr,"err = %g + %gi, ett = %g + %gi, eff = %g + %gi, ert = %g + %gi \n",
					creal(e_rr),cimag(e_rr), creal(e_tt),cimag(e_tt),
					creal(e_ff),cimag(e_ff), creal(e_rt),cimag(e_rt));
				  fprintf(stderr,"srr = %g + %gi, stt = %g + %gi, sff = %g + %gi, srt = %g + %gi \n",
					creal(s_rr),cimag(s_rr), creal(s_tt),cimag(s_tt),
					creal(s_ff),cimag(s_ff), creal(s_rt),cimag(s_rt));
				}
*/	
		/* Dissipated energy/vol based on strain */
		/* Imaginary component of stress*conj(strain) */
		work[e]	+= -(creal(s_rr) * cimag(e_rr) - creal(e_rr) * cimag(s_rr));
		work[e] += -(creal(s_tt) * cimag(e_tt) - creal(e_tt) * cimag(s_tt));
		work[e]	+= -(creal(s_ff) * cimag(e_ff) - creal(e_ff) * cimag(s_ff));
		work[e]	+= -(creal(s_rt) * cimag(e_rt) - creal(e_rt) * cimag(s_rt));
		work[e] /= 2.0;
		if (work[e] <= 0.0)
			work[e] == 0.0;
/*		if(e==1) fprintf(stderr,"work = %g",work[e]);*/

		/* Compute global dissipation */
		global_diss += 2.0*M_PI*work[e]*E->eco[e].area*pow(E->sphere.ro_dim,3)
											*E->tidal_visc[e];
		e_volume += 2.0*M_PI*E->eco[e].area;

		/* Volumetric Heating rate = E/V / dt */
		if(E->advection.timestep > 0.0)
			work[e] /= (dt);
		else
			work[e] = 0.0;
/*		if(e==1) fprintf(stderr,"--> %g",work[e]);*/

/*		fprintf(fp1,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
				creal(u_r),cimag(u_r),creal(u_th),cimag(u_th),
				creal(e_rr),cimag(e_rr),creal(e_tt),cimag(e_tt),
				creal(e_ff),cimag(e_ff),creal(e_rt),cimag(e_rt),
				creal(s_rr),cimag(s_rr),creal(s_tt),cimag(s_tt),
				creal(s_ff),cimag(s_ff),creal(s_rt),cimag(s_rt),
				work[e]);
  */
    
		/* Nondimensionalize Heating */
		work[e] *= E->sphere.ro_dim*E->sphere.ro_dim
						    /(E->data.density*E->data.Cp*E->data.DeltaT*E->data.therm_diff);

/*		if(e==1)  fprintf(stderr,"--> %g\n",work[e]);*/

  }

  /*fclose(fp1);*/
    
  /* Normalize to dissipation from conservation */
  for (e=1;e<=E->mesh.nel;e++)  {
    if ((E->control.despun == 0)  && (global_diss != 0.0))
			E->heating_despin[e] = work[e] * (Ediss/global_diss);
    else
			E->heating_despin[e] = 0.0;
/*		fprintf(E->fp,"e %d work %g heating %g\n",e,work[e],E->heating_despin[e]);*/
	/*	E->heating_despin[e] = 0.0;*/
  }
	fprintf(stderr,"Global Dissipation %g --> ",global_diss);
  global_diss *= 1.0/(E->data.density*E->data.Cp*E->data.DeltaT
																			*E->data.therm_diff*E->sphere.ro_dim);
	fprintf(stderr,"%g Volume %g Ediss %g\n\n",
				    global_diss,e_volume,Ediss);
	
    
  free((void *) log_c_Ri ); 
  free((void *) work ); 
}

/********************************************************
 * impact_heating                                       *
 *      function to simulate the deep mantle heating    *
 *      produced by a large impact.  Immediatley        *
 *      elevates the temperature within an isobaric     *
 *      core of radius R_ic by an amount dT.  dT decays *
 *      outside the core. DT is computed using the 			*
 *			shock heating model of Watters et al. (2009).  	*
 *			Calculating the approximation DT(P) as this is 	*
 *			more accurate than DT(u).                       *
 *                                                      *
 *			Now also calculating buoyancy due to isostatic	*
 *			adjustment of crater after material excavated.		*
 *																											*
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *      i       which impact in the list is this?       *
 ********************************************************/
void impact_heating(E,i)
        struct All_variables *E;
        int i;
{

	float xi,yi,zi;       /* Cartesian coords of hypocenter */
  float xn,yn,zn;       /* Cartesian coords of node */
  float delta;          /* Dist. btw. node and hypocenter */
  float dT,dT1;					/* Temperature increase */
	float Ed;							/* Energy density */

	int   e,j,k,m,n;
  float simp_comp;      /* Simple-to-complex transition */
  float D_at;           /* Transient crater diameter */
  float D_imp;          /* Diameter of impactor */
  float r_ic;           /* radius of isobaric core */
  float depth;          /* depth of center of isobaric core */
  float r_m;            /* radius of melt region */
  float V_m;            /* volume of melt region */
  float V_imp;          /* volume of impactor */
  float KE;             /* kinetic energy of impactor */
  float gamma;          /* efficiency of kinetic to thermal energy conversion */
	float H_at;						/* Transient crater depth at center */
	float *H_t;						/* Transient depth at all points in crater */
	float crater_sph_rad;	/* Radius of spherical "crater" */
	float slope_at;				/* slope of trans crater from horiz */
	float slope_atc;			/* Value of slope_at at rim */
	float rtheta;

	float latent;         /* specific internal energy of melting */
  float fsphere;        /* fraction of a sphere containing melt */
  float fm;             /* Vol. heated / Vol. IC */
  float A;              /* dT matching condition at r_m */
  float a1,b1,a2,b2,a,b,mvelo;  /* constants */
  float C,S,beta,fP,nP;         /* shock constants */
  float P_delta;                /* Shock pressure */
  float P_lith;         /* Lithostatic pressure */
  float u, uc;                  /* Particle velocity */
  float rcr;
	float F,Cmelt;							/* Melt fraction */
  int melt;             /* Flag for melted element */

  char output_file[255];
  FILE *fp1;

	float melting();
	float melting_by_node();
	float bulk_comp();
  void external_heating();
	void p_to_centres();
	void reallocate_markers();
	void get_markers_from_C();

  /*Set flag that impact is happening now */
  E->impacts.now = i;

  /* Crater data */

  simp_comp = 7.0e3;
  /*D_at = 0.87 * pow(simp_comp,0.115) *  pow(E->impacts.size[i]*1.0e3,0.885);*/ /* McKinnon & Schenk, 1985 */
  D_at = 0.98 * pow(simp_comp,0.079) *  pow(E->impacts.size[i]*1.0e3,0.92); /* Holsapple, 1993 */

	/* Trans. crater geometry.  This is sloppy for now */
	H_at = D_at/3.0;
	slope_atc = 1.176;
	crater_sph_rad = 1.624*H_at;

  D_imp = 0.69 * pow(D_at,1.27) * pow(E->data.grav_acc,0.28) 
                * pow(E->impacts.v,-0.56);
  V_imp = 4.0*M_PI/3.0 * pow(D_imp/2.0,3.0);

  r_ic = 0.5*D_imp * pow(10.0,a1) * pow(1.0e-3*E->impacts.v,b1);
  depth = 0.5*D_imp * pow(10.0,a2) * pow(1.0e-3*E->impacts.v,b2);

 /* Nondim  only need for Melosh*/
  depth *= 1.0/(E->sphere.ro_dim);
  r_ic *= 1.0/(E->sphere.ro_dim);
  r_m *= 1.0/(E->sphere.ro_dim);
  A = pow((r_ic/r_m),-4.4);

  /* Based of Melosh book, 1989 */
  latent = 3.4e6;
  gamma = 0.3;
  fm = 2.7;

/* Uncomment for Melosh heating
  r_ic = 0.72*D_imp;
  depth = r_ic;
*/
  V_m = 0.14 * E->impacts.v*E->impacts.v / latent;
  fsphere = 1.03 - 0.00195*V_m + 4.8e-4*V_m*V_m - 4.3e-6*V_m*V_m*V_m;
  V_m *= V_imp;
  r_m = pow((V_m * 3.0/(4.0*M_PI) * fsphere),(1.0/3.0));

  KE = 0.5*E->data.density*V_imp*E->impacts.v*E->impacts.v;
  KE *= gamma;
  KE -= latent*V_m*E->data.density;

  sprintf(output_file,"%s/impact.%d",E->control.data_file,i);
  fp1=fopen(output_file,"w");

  /* Get cartesian coords for impact hypocenter */
	/* In an axisymmetric case, an impact should really only go at the pole. 
		 Anywhere else and it's not circular */

	/* This part is only for 3D code.*/
	/*xi = (1.0-depth)*sin(E->impacts.th[i])*cos(E->impacts.f[i]);
  yi = (1.0-depth)*sin(E->impacts.th[i])*sin(E->impacts.f[i]);*/

	/* For axisymmetric, latitude is only surface direction and ought to only
	   be 0 or Pi */
  xi = (1.0-depth)*sin(E->impacts.th[i]);
	yi = 0.0;
  zi = (1.0-depth)*cos(E->impacts.th[i]);

	/* Don't bother to print phi or y */
  (void)fprintf(fp1,"%f %f %f %f %f %f %e %e\n",E->impacts.th[i],depth,xi,zi,
			E->impacts.size[i],r_ic,D_at,H_at);
  (void)fprintf(fp1,"%.6e %.6e %.6e %.6e\n",D_imp,V_imp,V_m,r_m);

  if(E->impacts.heat_from_file)
    external_heating(E,i,fp1);      /* Read heating from supplied files */
  else { /* compute impact heating following
          * Melosh/Pierazzo/Watters method */

	  /* Based off Watters et al. (2008) */

    a1 = -0.346;
    b1 = 0.211;
    a2 = -0.516;
    b2 = 0.361;
    a = -0.31;
    b = 1.15;
    mvelo = a + b*log10(1.0e-3*E->impacts.v);     /* Velocity decay exponent */

    C = 7.4e3;    /* intercept (m/s) */
    S = 1.25;     /* slope */
    beta = 0.5*C*C*E->data.density / S;
    uc = 0.5*E->impacts.v;
    nP = -1.84 + 2.61*log10(1.0e-3*E->impacts.v);

    (void)fprintf(fp1,"%.6e %.6e\n",nP,beta);

    /* Check each node */  
	  for(n=1;n<=E->mesh.nno;n++) {
      /* Get cartesian coords for node */
      /*xn = E->X[2][n]*sin(E->X[1][n])*cos(E->X[3][n]);
        yn = E->X[2][n]*sin(E->X[1][n])*sin(E->X[3][n]);*/
		  xn = E->X[2][n]*sin(E->X[1][n]);
      zn = E->X[2][n]*cos(E->X[1][n]);

      /*delta = sqrt((xn-xi)*(xn-xi) + (yn-yi)*(yn-yi) + (zn-zi)*(zn-zi));*/
      delta = sqrt((xn-xi)*(xn-xi) + (zn-zi)*(zn-zi));

      i = ((n-1) / E->mesh.noz) + 1;
      k = ((n-1) % E->mesh.noz) + 1;

      P_lith = E->data.density * E->data.grav_acc * (1.0-E->X[2][n]) 
              * E->sphere.ro_dim; /* lithostatic pressure */

		  /* Determine crater depth at all points */
		  if(k == E->mesh.noz){
			  rtheta = E->sphere.ro_dim * E->X[1][n]; /* dist. from impact */
			  if(rtheta < D_at/2.0) {
				  slope_at = asin(rtheta/crater_sph_rad);
				  E->impacts.H_t[i] = crater_sph_rad*(1.0 - cos(slope_at - slope_atc));
			  }
			  else {
				  E->impacts.H_t[i] = 0.0;
        }

			  fprintf(E->fpdebug,"i %d sl %f Ht %f\n",i,slope_at,E->impacts.H_t[i]);
		  }

		  if(delta < r_ic ) {	/* if the node is within the isobaric core */
        rcr = 1.0;				/* full effect of shock pressure*/
      } 
		  else {
      /*rcr = pow(r_ic/delta,1.87);*/
        rcr = pow(r_ic/delta,nP);
      }
      P_delta = E->data.density * (C + S*uc) * uc * rcr;
      /* P_delta = E->data.density * (C + S*uc*rcr) * uc*rcr; */
      /* For Foundering shock, subtract lithostatic pressure */
      P_delta -= P_lith;

      fP = (-P_delta/beta) / (1.0 - sqrt((2.0*P_delta/beta) +1) );
      dT = ( (P_delta/(2.0*E->data.density*S)) * (1 - 1.0/fP)
							- (C*C/(S*S)) * (fP - log(fP) - 1.0) )
							               / (E->data.Cp * E->data.ref_temperature);

		  /* dT is change in temperature, assuming no melting.*/
    }

	  if (E->control.melting) {
		  /*	F = melting_by_node(E,dT,n,k);
			E->Fm[n] = F; */
			(void)fprintf(fp1,"%f %f %f %f %e %e %.3f %.3f\n",E->X[1][n],
											E->X[2][n],delta,rcr,P_delta,dT,F,E->C[n]);
                      
      E->T[n] += dT;  /* We'll reduce this in a bit due to latent heat */
      /* By pre-adding dT here, we pass in dT=0 in call to melting.*/
		}
		else {
			/* Don't allow temperature to exceed the solidus */
			dT1 = dT;
			melt = 0;
//		(void)fprintf(stderr,"ifs %d %f\n",k,(E->T[n] + dT));
			if((E->T[n] + dT) > E->solidus[k]) {
				dT1 = max((E->solidus[k] - E->T[n]), 0.0);
				melt = 1;
//			(void)fprintf(stderr,"melt %d\n",melt);
			} 
			E->T[n] += dT1;
			(void)fprintf(fp1,"%f %f %f %f %e %e %e %d\n",E->X[1][n],E->X[2][n],
											delta,rcr,P_delta,dT,dT1,melt);
		}

//  fprintf(stderr,"i k n %d %d %d\n",i,k,n);
	} /* end internally generated impact heating */

	if (E->control.composition && E->control.melting) {
	  Cmelt = bulk_comp(E,5);

    /* Need T on element grid to do melting */
		p_to_centres(E,E->T,E->TE,E->mesh.levmax);

    /* Loop over elements to get melting by element */
		for(e=1;e<=E->mesh.nel;e++) {
			k = ((e-1) % E->mesh.elz) + 1;
			E->FmE[i] = melting(E,0,0,e,k);
		}

    /* Loop over nodes to reduce temperature */
		for(n=1;n<=E->mesh.nno;n++) {
			k = ((n-1) % E->mesh.noz) + 1;
			E->Fm[i] = melting_by_node(E,0,0,i,n);
		}

		/* Interpolate C,F from nodes to elements */
    /* Don't do this ! */
		/* p_to_centres(E,E->C,E->CE,E->mesh.levmax);
		 *  p_to_centres(E,E->Fm,E->FmE,E->mesh.levmax);
     */

		/* Update melt mass */
		E->Total.bulk_comp = 0.0;

		for(j=1;j<=E->mesh.esf;j++) { 
			E->slice.impact_melt[j] = 0.0;

			for(k=1;k<=E->mesh.elz;k++) {
				e = j*E->mesh.elz + k;
			//	E->slice.melt[j] += E->FmE[e]*E->eco[e].area/E->Total.vol;
			//	E->slice.new_melt[j] += E->FmE[e]*E->eco[e].area/E->Total.vol;
				E->slice.impact_melt[j] += E->FmE[e]*E->eco[e].area;
				E->Total.bulk_comp += E->CE[e] * E->eco[e].area;
			}

			E->slice.impact_melt[j] /= E->Total.vol;
			E->slice.melt[j] +=	E->slice.impact_melt[j];
			E->slice.new_melt[j] +=	E->slice.impact_melt[j];
		}
		
		E->Total.bulk_comp /= E->Total.vol;

	  Cmelt = bulk_comp(E,6);
	}
	(void)fclose(fp1);

	return;

}


/********************************************************
 * external_heating                                     *
 *      function to read in the mantle heating          *
 *      produced by a large impact from an external     *
 *      file. Interpolates from regular grid onto       *
 *      Citcom grid.  Heating from external file is     *
 *      axisymmetric.                                   *
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *      i       which impact in the list is this?       *
 *      fp1     impact heating output file              *
 ********************************************************/
void external_heating(E,i,fp1)
        struct All_variables *E;
        int i;
        FILE *fp1;
{
  float dT,dT1;
  float *theta_in, *r_in;       /* Regular grid coords */
  float dr_in, dth_in;          /* Regular Grid intervals */
  float *dT_in;                 /* T increase on reg. grid */
  float dist;                   /* distance btw CitcomS grid pt and impact axis */
  float A_ll,A_lr,A_ul,A_ur,A_nodes;            /* Weighting functions for interpolation */
  float temp1,temp2,temp3;

  int num_theta, num_r;         /* Number of input grid points in latitude and radius */
  int nodes_in;
  int j,k,m,n;          /* Counters */
  int left, right, lower, upper;
  int ul,ur,ll,lr;              /* Interpolants; coords of nearest neighbors */
  int melt;                     /* Flag for melting */
  int rnode;

  char input_file[255];
  char input_s[200];
  FILE *fp0;

  /*
   * First, read in heating on regular grid.
   */
  fprintf(stderr,">>> %s\n",E->impacts.heating_file);
  sprintf(input_file,"%s/%s",E->control.data_file,E->impacts.heating_file);
  fprintf(stderr,"%s\n",input_file);
  fp0=fopen(input_file,"r");

  /* Should stick this in input file */
  num_theta = 91;
  num_r = 41;
  nodes_in = num_theta*num_r;

  dth_in = M_PI / (float)(num_theta-1);
  dr_in = (E->sphere.ro - E->sphere.ri) / (float)(num_r-1);

  theta_in = (float *)malloc((num_theta+1)*sizeof(float));
  r_in = (float *)malloc((num_r+1)*sizeof(float));
  dT_in = (float *)malloc((nodes_in+1)*sizeof(float));

  for(j=1;j<=num_theta;j++)
    theta_in[j] = dth_in*(j-1);

  for(k=1;k<=num_r;k++)
    r_in[k] = E->sphere.ri + dr_in*(k-1);

  r_in[1] -= 1.0e-7;    /* Make sure lowest point is not above CMB */
  r_in[num_r] += 1.0e-7;        /* Make sure highest point is not below surface */

  for(j=1;j<=num_theta;j++) { 
    for(k=1;k<=num_r;k++) {
      n = (j-1)*num_r + k;
      fgets(input_s,100,fp0);
      sscanf(input_s,"%f ",&(dT_in[n]));
      if (dT_in[n] >= 5.0e4)
        dT_in[n] = 0.0;
      dT_in[n] *= 0.9;  // Reduce for testing
    fprintf(E->fp,"n %d theta %f r %f dTin %f\n",n,theta_in[j],r_in[k],dT_in[n]);
      }
   }

  (void)fclose(fp0);

  /*
   * Interpolate dT onto CitcomS grid 
   */
  
  /* Visit each CitcomS node */
  for(n=1;n<=E->mesh.nno;n++) {

      rnode = (n-1) % E->mesh.noz + 1; /* Radial position of each node */
      //fprintf(E->fp,"n %d rnode %d\n",n,rnode);

      /* Get lat distance btw node and impact axis */
      /* For a north polar impact, dist = E->X[1][n]; */
      dist = E->X[1][n] - E->impacts.th[i];
//      fprintf(stderr,"dist %f\n",dist);

     /*
      * Scan the regular grid to find nearest neighbors
      * Probably there is a more efficient way to do this,
      * but it's only done once, so who cares.
      */
      left = right = lower = upper = 0;
      for (j=1; j<num_theta; j++) {
      /* If the citcom elem is between two adjacent regular nodes horiz. */
        if((theta_in[j] <= dist) && (theta_in[j+1] >= dist)) {
          left = j;
          right = j+1;
        }
      }
      for (k=1; k<num_r; k++) {
      /* and vertically */
        /*if (n==65) {
          fprintf(E->fp_out,"k %d %g %g %g\n",k,r_in[k],E->X[2][n],r_in[k+1]);
          fprintf(E->fp_out,"k %d %g\n",k,r_in[k]-E->X[2][n]);
        }*/
        if((r_in[k] <= E->X[2][n]) && (r_in[k+1] >= E->X[2][n])) {
          lower = k;
          upper = k+1;
        }
      }

      ll = (left-1)*num_r + lower;
      ul = (left-1)*num_r + lower + 1;
      lr = left*num_r + lower;
      ur = left*num_r + lower + 1;


      /*
       * Average the heating values at the four regular elements adjacent to
       * the Citcom element in question, weighted by their proximity.
       */

      A_ll = (theta_in[left] - dist) * (r_in[lower] - E->X[2][n]);
      A_ul = (theta_in[left] - dist) * -(r_in[upper] - E->X[2][n]);
      A_lr = -(theta_in[right] - dist) * (r_in[lower] - E->X[2][n]);
      A_ur = -(theta_in[right] - dist) * -(r_in[upper] - E->X[2][n]);

      /* Total Area of quadrangle formed by adjacent elements */
      A_nodes = A_ll + A_ul + A_lr + A_ur;

      /* Weighting is proportional to the area of the opposite regular node */
      /* Larger area means node is farther away.  Lever rule. */
      dT = (dT_in[ll]*A_ur + dT_in[lr]*A_ul + dT_in[ul]*A_lr + dT_in[ur]*A_ll) / A_nodes;
      /* dT is change in temperature, assuming no melting.*/
      dT /= E->data.ref_temperature;  /* Nondimensionalize */

      if (E->control.melting) {
         (void)fprintf(fp1,"%f %f %e %.3f %.3f\n",E->X[1][n],E->X[2][n],dist,
                        dT,E->C[n]);
         E->T[n] += dT;  /* We'll reduce this in a bit due to latent heat */
         /* By pre-adding dT here, we pass in dT=0 in call to melting.*/
      }
      else {
        /* Don't allow temperature to exceed the solidus */
        dT1 = dT;
        melt = 0;
        if((E->T[n] + dT) > E->solidus[rnode]) {
          dT1 = max((E->solidus[rnode] - E->T[n]), 0.0);
          melt = 1;
        (void)fprintf(stderr,"ifs %d %f %f\n",rnode,(E->T[n] + dT),E->solidus[rnode]);
                (void)fprintf(stderr,"melt %d\n",melt);
        }
        E->T[n] += dT1;
        (void)fprintf(fp1,"%f %f %f %e %e %d\n",E->X[1][n],E->X[2][n],dist,dT,
               dT1,melt);
      } 
    }

  return;

}

/********************************************************
 * bulk_comp			                                      *
 *      function to compute the bulk composition of the *
 *			mantle.																					*
 *																											*
 * Parameters                                           *
 *      E       All_variables                           *
 *		  i				Index for debugging											*
 *																											*
 * Returns																							*
 *			Globally averaged composition										*
 *																											*
 ********************************************************/

float bulk_comp(E,i) 
  struct All_variables *E;
	int i;
{
	float comp; /* bulk composition */
	int		e;		/* element counter */
	void p_to_centres();

//	p_to_centres(E,E->C,E->CE_temp,E->mesh.levmax);

  comp = 0.0;

	for(e=1;e<=E->mesh.nel;e++){
		//comp += E->CE_temp[e] * E->eco[e].area;
		comp += E->CE[e] * E->eco[e].area;
		if (E->advection.timesteps <= 25  && e == 9216) {
			fprintf(E->fpdebug,"e %d %g %g\n",e,E->CE[e],comp);
	    fprintf(E->fpdebug,"> Cn = %g %g %g %g\n",
              E->C[9407],E->C[9408],E->C[9456],E->C[9457]);
      }
	}
			
	comp /= E->Total.vol;
	fprintf(E->fpdebug,"> %d bulk comp = %g %g \n",i,comp,E->Total.bulk_comp);
//	fprintf(E->fpdebug,"> tracer = %g\n",E->C12f[1806336]);

	return(comp);

}
