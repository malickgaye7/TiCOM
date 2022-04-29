/*  Here are the routines which process the results of each buoyancy solution, and call
    any relevant output routines. Much of the information has probably been output along
    with the velocity field. (So the velocity vectors and other data are fully in sync).
    However, heat fluxes and temperature averages are calculated here (even when they
    get output the next time around the velocity solver);
    */

#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif
//#include <malloc.h>
#include <sys/types.h>
//#include <stdlib.h> /* for "system" command */
//#include <malloc.h>

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#include <stdio.h>

#ifndef __ELEMENT_DEFINITIONS_H__
#define __ELEMENT_DEFINITIONS_H__
#include "element_definitions.h"
#endif

#ifndef __GLOBAL_DEFS_H__
#define __GLOBAL_DEFS_H__
#include "global_defs.h"
#endif

void process_temp_field(E,ii)
 struct All_variables *E;
    int ii;
{ 
    void heat_flux();
    void output_temp();
    float freezemelt();

    float Fs, Fb;     /* Global surface and bottom heat fluxes */
    float deltah;     /* Change in shell thickness this timestep */
    int   i,j;        /* Node indices */

        /* Compute heat flux */
//    if ( ((ii % E->control.record_every) == 0))    {
        heat_flux(E);
        
        Fs=0.0;
        for (i=1;i<=E->mesh.nox;i++)  {
          j=i*E->mesh.noz;
          if (i>1) {
            Fs += (E->slice.shflux[i]+E->slice.shflux[i-1])*0.5*
                sin(0.5*(E->X[1][j]+E->X[1][j-E->mesh.noz]))*
                (E->X[1][j]-E->X[1][j-E->mesh.noz]);
          }
        }
        Fs = Fs/(1.0-cos(E->X[1][E->mesh.nno]));
        
        Fb=0.0;
        for (i=1;i<=E->mesh.nox;i++)  {
          j=i*E->mesh.noz;
          if (i>1) {
            Fb += (E->slice.bhflux[i]+E->slice.bhflux[i-1])*0.5*
                sin(0.5*(E->X[1][j]+E->X[1][j-E->mesh.noz]))*
                (E->X[1][j]-E->X[1][j-E->mesh.noz]);
          }
        }
        Fb = Fb/(1.0-cos(E->X[1][E->mesh.nno]));


/*
      output_temp(E,ii);
*/
//      }

      /* Compute thickening or thinning of ice shell */
      if(E->control.freezing == 1) {
        deltah = freezemelt(E,Fb);
        //fprintf(stderr,"%e\n",deltah);
        E->sphere.deltarb += deltah; /* Add to running total this regrid */
        E->monitor.deltah += deltah; /* Add to running total */

        fprintf(stderr,"%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n",
                    E->monitor.elapsed_time,Fs,Fb,deltah,E->sphere.deltarb,
                    E->monitor.deltah);
        fprintf(E->fpq,"%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n",
                    E->monitor.elapsed_time,Fs,Fb,deltah,E->sphere.deltarb,
                    E->monitor.deltah);
        fflush(E->fpq);
        }

    return;
}
/* ===================
    Surface heat flux  
   =================== */

void heat_flux(E)
    struct All_variables *E;
{
    int ee,e,i,j,node,lnode;
    float *mass,*flux,*SU,*RU;
    float diff,T1[9],VZ[9],u[9],T[9],dTdz[9],area,uT,adv;

    double xk[3][5];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    void get_global_shape_fn();
    void return_horiz_ave();

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->mesh.nno;
    const int lev = E->mesh.levmax;

    flux = (float*) malloc((nno+1)*sizeof(float));
    return_horiz_ave(E,E->T,E->Have.T);

    for(i=1;i<=nno;i++) {
      flux[i] = 0.0;
      E->heatflux[i] = 0.0;
    }
    for(e=1;e<=E->mesh.nel;e++) {
      ee = (e-1)%E->mesh.elz+1;
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,xk,2,E->mesh.levmax);
      for(j=1;j<=ends;j++) {
        VZ[j] = E->V[2][E->ien[e].node[j]];
      }
      for(i=1;i<=ppts;i++)   {
        u[i] = 0.0;
        T[i] = 0.0;
        T1[i] = 0.0;
        dTdz[i] = 0.0;
        for(j=1;j<=ends;j++)  {
          lnode = (E->ien[e].node[j]-1)%E->mesh.noz+1;
          u[i] += VZ[j]*E->N.ppt[GNPINDEX(j,i)];
          T[i] += E->T[E->ien[e].node[j]]*E->N.ppt[GNPINDEX(j,i)];
          T1[i] += (E->T[E->ien[e].node[j]]-E->Have.T[lnode])*E->N.ppt[GNPINDEX(j,i)];
          dTdz[i] += -E->T[E->ien[e].node[j]]*GNx.ppt[GNPXINDEX(1,j,i)];
        }
      }
      uT = 0.0;
      area = 0.0;
      adv = 0.0;
      diff = (E->diffusivity[ee]+E->diffusivity[ee+1])*0.5;
      for(i=1;i<=ppts;i++)   {
        uT += u[i]*T[i]*dOmega.ppt[i] + diff*dTdz[i]*dOmega.ppt[i];
        adv += u[i]*T1[i]*dOmega.ppt[i];
        area += dOmega.ppt[i];
      }
      uT /= area;
      adv /= area;
      for(j=1;j<=ends;j++)  {
        flux[E->ien[e].node[j]] += uT*E->TWW[E->mesh.levmax][e].node[j];
        E->heatflux[E->ien[e].node[j]] += adv*E->TWW[E->mesh.levmax][e].node[j];
      }
    }           /* end of e */

    for(i=1;i<=nno;i++)   {
      flux[i] = flux[i]*E->Mass[i];
      E->heatflux[i] = E->heatflux[i]*E->Mass[i];
    }

    for(i=1;i<=E->mesh.nsf;i++)   {
      E->slice.shflux[i] = 2*flux[E->surf_node[i]]
                           - flux[E->surf_node[i]-1];

      E->slice.bhflux[i] = 2*flux[E->surf_node[i]-E->mesh.noz+1]
                           - flux[E->surf_node[i]-E->mesh.noz+2];
    }

    return_horiz_ave(E,flux,E->Have.f);
    return_horiz_ave(E,E->heatflux,E->Have.F);


   free(flux);

  return;  
  }
  