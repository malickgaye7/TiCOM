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
#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif
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

void general_stokes_solver(E)
     struct All_variables *E;

{
    void construct_stiffness_B_matrix();
    void velocities_conform_bcs();
    void assemble_forces();
    float fvdot();
    float vnorm_nonnewt();
    void get_system_viscosity();
    void solve_constrained_flow_iterative();

    float vmag;
    double Udot_mag = 0.0;
    double dUdot_mag = 0.0;
    double CPU_time0(),time;
    double  kineticE();
    int count,i,j,k;

    static float *oldU,*delta_U;
    static int visits=0;

    const int nno = E->mesh.nno;
    const int nel = E->mesh.nel;
    const int nnov = E->mesh.nnov;
    const int neq = E->mesh.neq;
    const int vpts = vpoints[E->mesh.nsd];
    const int dims = E->mesh.nsd;
    const int addi_dof = additional_dof[dims];

    if(visits==0) {
	    oldU = (float *)malloc((neq+2)*sizeof(float));
      delta_U = (float *)malloc((neq+2)*sizeof(float));
	    for(i=0;i<=neq;i++) 
	      oldU[i]=0.0;
      visits++;
    }

     
    /* FIRST store the old velocity field */

    E->monitor.elapsed_time_vsoln1 =  E->monitor.elapsed_time_vsoln;
    E->monitor.elapsed_time_vsoln = E->monitor.elapsed_time;

    time=CPU_time0();

    velocities_conform_bcs(E,E->U);

    assemble_forces(E,0); 

    count=1;
    
    do  {

      if(E->viscosity.update_allowed)
          get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);

      construct_stiffness_B_matrix(E);

      solve_constrained_flow_iterative(E);	

      Udot_mag = kineticE(E,E->U,E->mesh.levmax);
      fprintf(E->fp,"%g %.6e \n",E->monitor.elapsed_time,Udot_mag);
      fflush(E->fp);

      //fprintf(stderr,"kinetic energy= %.7e at time= %g for step %d\n",Udot_mag,E->monitor.elapsed_time,E->monitor.solution_cycles);
      if (E->control.dimensionalize) {
        double dt = E->monitor.elapsed_time * E->monitor.time_scale;
        fprintf(stderr, "%5d\t%20g\t\t %4.7e\n", E->monitor.solution_cycles,dt,Udot_mag);
      } else { fprintf(stderr, "%5d\t%20g\t\t %4.7e\n", E->monitor.solution_cycles,E->monitor.elapsed_time,Udot_mag); }


      if (  E->viscosity.SDEPV  )   {
        for (i=0;i<neq;i++) {
          delta_U[i] = E->U[i] - oldU[i]; 
          oldU[i] = E->U[i];
          }
        Udot_mag  = sqrt(fvdot(E,oldU,oldU,E->mesh.levmax));
        dUdot_mag = vnorm_nonnewt(E,delta_U,oldU,E->mesh.levmax); 

          fprintf(stderr,"Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d\n",dUdot_mag,Udot_mag,count);
          fprintf(E->fp,"Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d\n",dUdot_mag,Udot_mag,count);
          fflush(E->fp); 
        count++;
        }         /* end for SDEPV   */

      } while((count < 50) && (dUdot_mag>E->viscosity.sdepv_misfit) && E->viscosity.SDEPV);

          fflush(E->fp); 

  return;
}
