 /*CITCOM: A finite element convection program written at Caltech 1992 */
 /*Aims to include an iterative matrix solver based on Multigrid techniques */
 /*To do this requires the use of a mixed method and a conjugate-gradient */
 /*approach to determining the */


// The 2D axisymmetric spherical version of citcom was implemented
// by Shijie Zhong in 1999. The particle method for composition was
// added in 2002 by SZ.
//
// A reference for this code is 
//  J. H. Roberts, and S. Zhong, Plume-induced topography and geoid anomalies
//   and their implications for the Tharsis Rise on Mars, JGR, 109, Art. No. E03009,
//   doi:10.1029/2003JE002226, 2004.
//
// If you use this code, you may reference this paper and also L. Moresi
// 's original paper on Cartesian citcom.
//
// The particle method works significantly better than the field method
// in my view. Although there is a version of field method implemented
// in this code, I do not recommend to use it.
//                          a note by SZ


#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#ifndef STRING_H
#define STRING_H
#include <string.h>
#endif
//#include <malloc.h>
#include <sys/types.h>

#ifndef __ELEMENT_DEFINITIONS_H__
#define __ELEMENT_DEFINITIONS_H__
#include "element_definitions.h"
#endif

#ifndef __GLOBAL_DEFS_H__
#define __GLOBAL_DEFS_H__
#include "global_defs.h"
#endif

extern int Emergency_stop;

main(argc,argv)
     int argc;
     char **argv;
     
{	/* Functions called by main*/
  void general_stokes_solver();
  void read_instructions();
  void solve_constrained_flow();
  void solve_derived_velocities();
  void process_temp_field(); 
  void process_heating();
	void reallocate_markers();
	void get_markers_from_C();
	
	float bulk_comp();
  float dot();

  int k, *temp;
  double CPU_time0(),time,initial_time,start_time;
	float C;
 
  struct All_variables E;

 E.monitor.solution_cycles=0;

  start_time = time = CPU_time0();
 
  read_instructions(&E,argc,argv);

  E.control.keep_going=1;

     fprintf(stderr,"Input parameters taken from file '%s'\n\n",argv[1]);
     fprintf(stdout,"Initialization complete after %g seconds\n\n",CPU_time0()-time); fflush(E.fp);
     initial_time = CPU_time0()-time;
     E.monitor.cpu_time_on_vp_it = CPU_time0();

  if (E.control.dimensionalize) { fprintf(stderr, "Units -> Time: Years\n\t Energy: Nondimensionalized\n\n"); }
  else { fprintf(stderr, "Units -> Nondimensionalized\n\n"); }
  fprintf(stderr, " \e[4mStep\t\tTime\t\t\tKinetic Energy\e[m\n");
  general_stokes_solver(&E);
  process_temp_field(&E,E.monitor.solution_cycles); 

  void process_new_velocity(struct All_variables *, int);
  process_new_velocity(&E,E.monitor.solution_cycles);

  if (E.control.stokes)  {
     E.control.keep_going=0;
     E.monitor.solution_cycles++; 
     }

  while ( E.control.keep_going   &&  (Emergency_stop == 0) )   {
		
			fprintf(E.fpdebug,"\nStep: %d\n",E.advection.timesteps);
      E.monitor.solution_cycles++; 
      if(E.monitor.solution_cycles>E.control.print_convergence)
         E.control.print_convergence=1;

      process_heating(&E);

			C = bulk_comp(&E,0);
      (E.next_buoyancy_field)(&E);

      process_temp_field(&E,E.monitor.solution_cycles); 
   
      if (E.control.freezing && ((E.sphere.deltarb > 0.5*(E.eco[1].size[2])) || E.control.keep_going==0) ) {
        void adjust_model_domain(struct All_variables *);
        adjust_model_domain(&E);
      }

      E.monitor.elapsed_time += E.advection.timestep;
      E.impacts.now = -1; /* Reset impact flag */

      //fprintf(stderr, " \e[4mStep\t\tTime\t\t\tKinetic Energy\e[m\n");
      general_stokes_solver(&E);

			if (E.control.composition && strcmp(E.control.comp_adv_method,"particle")==0) {
      	(E.next_buoyancy_field)(&E);      /* correct with R-G */
      }

      if (E.control.composition && strcmp(E.control.comp_adv_method,"particle")==0) {

			  /* Assign tracers based on composition */
			//	reallocate_markers(&E);

				get_markers_from_C(&E,E.CE,E.CElement);
				
				(void)fprintf(stderr,"Markers curr: %d %d %d\n",
								E.advection.marker_type[0],E.advection.marker_type[1],
								E.advection.markers);
				(void)fprintf(E.fpdebug,"Markers curr: %d %d\n",
								E.advection.marker_type[0],E.advection.marker_type[1]);
		
				 }

      process_new_velocity(&E,E.monitor.solution_cycles);

        fprintf(E.fp,"CPU total = %g & CPU = %g for step %d time = %.4e dt = %.4e  maxT = %.4e sub_iteration%d\n",CPU_time0()-start_time,CPU_time0()-time,E.monitor.solution_cycles,E.monitor.elapsed_time,E.advection.timestep,E.monitor.T_interior,E.advection.last_sub_iterations);
        fflush(E.fp);

        time = CPU_time0();

      }
  
     E.monitor.cpu_time_on_vp_it=CPU_time0()-E.monitor.cpu_time_on_vp_it;
     fprintf(E.fp,"Initialization overhead = %f\n",initial_time);
     fprintf(E.fp,"Average cpu time taken for velocity step = %f\n",
	 E.monitor.cpu_time_on_vp_it/((float)(E.monitor.solution_cycles)));
     fprintf(stdout,"Average cpu time taken for velocity step = %f\n",
	 E.monitor.cpu_time_on_vp_it/((float)(E.monitor.solution_cycles)));

  fclose(E.fp);
  // DEBUG
  //deallocate_common_vars(E);
  int i,j;
  //for(i=E.mesh.levmin;i<=E.mesh.levmax;i++)  {
    free(E.surf_node);
    free(E.surf_element);
    free(E.sien);
  //}
    //
  for(i=E.mesh.levmin;i<=E.mesh.levmax;i++)   {
    free(E.NEI[i].element);
    free(E.NEI[i].lnode);
    free(E.NEI[i].nels);

    free(E.TWW[i]);
    free(E.NODE[i]);
    free(E.VI[i]);
    free(E.TW[i]);

    free(E.EVI[i]);

    free(E.BPI[i]);
    free(E.elt_del[i]);
    free(E.LMD[i]);
    free(E.EL[i]);
    free(E.IEN[i]);
    free(E.ID[i]);
    free(E.ECO[i]);

    free(E.MASS[i]);

    for(j=1;j<=E.mesh.nsd;j++)   {
      free(E.Interp[i][j]);
      free(E.XX[i][j]);
    }
  }
  
  free(E.XP[2]);
  free(E.XP[1]);
  free(E.liquidus);
  free(E.lherzliq);
  free(E.solidus);
  free(E.expansivity);
  free(E.diffusivity);
  free(E.mat);
  //free(E->Total.melt_prod);
  free(E.slice.impact_melt);
  free(E.slice.new_melt);
  free(E.slice.melt);
  free(E.slice.vxsurf[2]);
  free(E.slice.vxsurf[1]);
  free(E.slice.cen_mflux);
  free(E.slice.bhflux);
  free(E.slice.shflux);
  free(E.slice.vlinek);
  free(E.slice.vline);
  free(E.slice.tpgb);
  free(E.slice.tpg);
  free(E.stress);
  free(E.Have.Tprev);
  free(E.Have.Tadi);
  free(E.Have.F);
  free(E.Have.f);
  free(E.Have.vrms);
  free(E.Have.Rho);
  free(E.Have.r_Ri);
  free(E.Have.i_Ri);
  free(E.Have.Vi);
  free(E.Have.T);
  for(i=1;i<=E.mesh.nsd;i++)  {
    free(E.TB[i]);
    free(E.CB[i]);
  }
  free(E.Fas410_b);
  free(E.Fas410);
  free(E.Fas670_b);
  free(E.Fas670);
  free(E.c_ERi);
  free(E.c_Ri);
  free(E.heating_despin);
  free(E.heating_shear);
  free(E.tidal_visc);
  free(E.heating_tidal);
  free(E.heating_adi);
  free(E.heating_latent);
  free(E.heating_visc);
  free(E.heatflux);
  free(E.edot);
  free(E.NP);
  free(E.buoyancy);
  free(E.TE);
  free(E.T);
  free(E.F);
  free(E.U);
  free(E.FmE);
  free(E.Fm);
  free(E.CE_temp);
  free(E.CE_prev);
  free(E.CE);
  free(E.C_prev);
  free(E.C);

  int lev;
  for(lev=E.mesh.levmax;lev>=E.mesh.levmin;lev--) {
    free(E.Node_map[lev]);
    free(E.Eqn_k[lev]);
  }

  free(E.S);
  free(E.P);
  free(E.Tdot);
  for(i=1;i<=E.mesh.nsd;i++)  {
      free(E.V[i]);
      free(E.VB[i]);
  }

  for(i=E.mesh.levmin;i<=E.mesh.levmax;i++)  {
    free(E.BI[i]); 
    free(E.EQN[i]); 
  }
  for (i=E.mesh.levmin;i<=E.mesh.levmax;i++) {
    if(!E.control.NMULTIGRID && !E.control.NASSEMBLE)  {
       free(E.elt_k[i]);
    }
    else if (E.control.NMULTIGRID || E.control.NASSEMBLE) {
      free(E.Node_eqn[i]);
      free(E.Node_k_id[i]);
    }
  }
  if (E.control.composition) { free(E.Cdot); }
  if (!(strcmp(E.control.comp_adv_method,"field")==0)) {
    for(i=1;i<=E.mesh.nsd;i++)   {
      free(E.VO[i]);
      free(E.XMC[i]);
      free(E.XMCpred[i]);
    }
    free(E.C12);
    free(E.C12f);
    free(E.CElement);
  }
  return 0;  
}
