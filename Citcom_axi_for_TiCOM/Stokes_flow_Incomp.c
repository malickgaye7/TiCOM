/*   Functions which solve for the velocity and pressure fields using Uzawa-type iteration loop.  */

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

extern int Emergency_stop;

/* Master loop for pressure and (hence) velocity field */


void solve_constrained_flow_iterative(E)
     struct All_variables *E;

{ 
    double *D1;
    double *u;
    double *R,*Bp;
    double residual_ddash;
    double vmag;
    double vdot(),pdot();

    static int been_here = 0;
   
    float solve_Ahat_p_fhat();
    void v_from_vector();
    void dp_to_nodes();   
   
    int steps,cycles;
    int i,j,k,doff,vel_cycles_previous,vel_calls_previous;
  
    double time,CPU_time0();
 
    const int npno = E->mesh.npno;
    const int gnpno = E->mesh.npno;
    const int nno = E->mesh.nno;
    const int dims = E->mesh.nsd;
    const int neq = E->mesh.neq;
    const int gneq = E->mesh.neq;

    time=CPU_time0();
    
    cycles=E->control.p_iterations;
   
    /* Solve for velocity and pressure, correct for bc's */

    residual_ddash=solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

    been_here=1;

    v_from_vector(E,E->V,E->U); 
    dp_to_nodes(E,E->P,E->NP,E->mesh.levmax);   

  return; 
}



/*  ==========================================================================  */

float solve_Ahat_p_fhat(E,V,P,F,imp,steps_max)

     struct All_variables *E;
     double *V,*P,*F;
     double imp;
     int *steps_max;
    
{ 
  int i,j,k,ii,count,convergent,valid,problems,lev,lev_low,npno,neq,steps;
  int gnpno,gneq;
  
  static int been_here = 0;
  double *p1,*r1,*u;
  double *r0,*r2,*z0,*z1,*s1,*s2,*Ah,*u1;
  double *shuffle, *R;
  double alpha,delta,s2dotAhat,r0dotr0,r1dotz1;
  double residual, initial_residual, last_residual,v_res;
  double dpressure,dvelocity;
  
  double vdot(),pdot();

  float CPU_time();
  double time0,time,CPU_time0();

  void assemble_div_u();
  void assemble_del2_u();
  void assemble_grad_p();
  void strip_bcs_from_residual();
  int  solve_del2_u();
   
  const int dims=E->mesh.nsd;
  const int n=loc_mat_size[E->mesh.nsd];

  npno=E->mesh.npno;
  neq=E->mesh.neq;

  gnpno=E->mesh.npno;
  gneq=E->mesh.neq;

  r0 = (double *)malloc((npno+1)*sizeof(double));
  r1 = (double *)malloc((npno+1)*sizeof(double));
  r2 = (double *)malloc((npno+1)*sizeof(double));
  z0 = (double *)malloc((npno+1)*sizeof(double));
  z1 = (double *)malloc((npno+1)*sizeof(double));
  s1 = (double *)malloc((npno+1)*sizeof(double));
  s2 = (double *)malloc((npno+1)*sizeof(double));
  p1 = (double *)malloc((npno+1)*sizeof(double));
  Ah = (double *)malloc((neq+1)*sizeof(double));
  u1 = (double *)malloc((neq+1)*sizeof(double));

  problems=0;
  time0=time=CPU_time0();

  been_here ++;
 
  /* calculate the velocity residual, note there are tricks involved here */

  lev=E->mesh.levmax;

  assemble_grad_p(E,P,Ah,lev);
  assemble_del2_u(E,V,u1,lev,1);

  for(i=0;i<neq;i++) 
      Ah[i] = F[i] - Ah[i] - u1[i]; 

  v_res=sqrt(vdot(E,F,F,lev)/gneq);

  //fprintf(stderr,"initial residue of momentum equation F %.8e %d\n",v_res,gneq);
  fprintf(stdout,"initial residue of momentum equation F %.8e %d\n",v_res,gneq); //DEBUG

  strip_bcs_from_residual(E,Ah,lev);

  valid=solve_del2_u(E,u1,Ah,imp*v_res,E->mesh.levmax);
  strip_bcs_from_residual(E,u1,lev);

  if(!valid) problems++; 
  
  for(i=0;i<neq;i++)  {
      V[i] += u1[i];
     }

  assemble_div_u(E,V,r1,lev);

  residual = initial_residual = sqrt(pdot(E,r1,r1,lev)/gnpno);

  E->monitor.vdotv = sqrt(vdot(E,V,V,lev)/gneq);

  E->monitor.incompressibility = residual/E->monitor.vdotv;
         
   for(i=1;i<=npno;i++)
        p1[i] = 0.0;
     
   count = 0;
   convergent=0;

   if (E->control.print_convergence)  {
         fprintf(E->fp,"AhatP (%03d) after %g seconds with div/v=%.3e for step %d\n",count,CPU_time0()-time0,E->monitor.incompressibility,E->monitor.solution_cycles); /**/
         fflush(E->fp);
         }         

   dpressure = 1.0;
   dvelocity = 1.0;

  while( count==0 || ((count < *steps_max) && (dvelocity >= imp || dpressure >=imp) && E->monitor.incompressibility >= E->control.tole_comp) )  { 
/*   while( count==0 || ((count < *steps_max) && (dvelocity >= imp || dpressure >=imp) ) )  {*/

     for(j=1;j<=npno;j++)
       z1[j] = E->BPI[lev][j]*r1[j];
     
     r1dotz1 = pdot(E,r1,z1,lev);

     if ((count == 0))
       for(j=1;j<=npno;j++)
            s2[j] = z1[j];
     else {
       r0dotr0=pdot(E,r0,z0,lev);
       assert(r0dotr0 != 0.0  /* Division by zero in head of incompressibility iteration */);
       delta = r1dotz1/r0dotr0;
       for(j=1;j<=npno;j++)
            s2[j] = z1[j] + delta * s1[j];
       }
      
     assemble_grad_p(E,s2,Ah,lev); 

     valid=solve_del2_u(E,u1,Ah,imp*v_res,lev);  
     strip_bcs_from_residual(E,u1,lev);

     if(!valid) problems++;
      
     assemble_div_u(E,u1,Ah,lev);

     s2dotAhat=pdot(E,s2,Ah,lev);
      
	                 /* alpha defined this way is the same as R&W */
     alpha = r1dotz1/s2dotAhat; 
     
     for(j=1;j<=npno;j++)   {
       r2[j] = r1[j] - alpha * Ah[j];
       p1[j] += alpha * s2[j];
       }
     
     for(j=0;j<neq;j++)
       V[j] -= alpha * u1[j];
      
     assemble_div_u(E,V,Ah,lev);
     E->monitor.vdotv = vdot(E,V,V,E->mesh.levmax);
     E->monitor.incompressibility = sqrt((gneq/gnpno)*(1.0e-32+pdot(E,Ah,Ah,lev)/(1.0e-32+E->monitor.vdotv)));

     dpressure=alpha*sqrt(pdot(E,s2,s2,lev)/(1.0e-32+pdot(E,p1,p1,lev)));
     dvelocity=alpha*sqrt(vdot(E,u1,u1,lev)/(1.0e-32+E->monitor.vdotv));
         
     count++;
     if (E->control.print_convergence && count%5==0)  {
       fprintf(E->fp,"AhatP (%03d) after %g seconds with div/v=%.3e for step %d dv=%g dp=%g\n",count,CPU_time0()-time0,E->monitor.incompressibility,E->monitor.solution_cycles,dvelocity,dpressure); /**/
       fflush(E->fp);
       }         

     shuffle=s1;s1=s2;s2=shuffle;
     shuffle=r0;r0=r1;r1=r2;r2=shuffle;
     shuffle=z0;z0=z1;z1=shuffle;

     }       /* end loop for conjugate gradient   */

    if(problems) {
      fprintf(E->fp,"Convergence of velocity solver may affect continuity\n");
      fprintf(E->fp,"Consider running with the `see_convergence=on' option\n");
      fprintf(E->fp,"To evaluate the performance of the current relaxation parameters\n");
      fflush(E->fp);
      }

  for(j=1;j<=npno;j++) 
      P[j] += p1[j];

    
  free((void *) r0);
  free((void *) r1);       
  free((void *) r2);
  free((void *) z0);
  free((void *) z1);
  free((void *) s1);
  free((void *) s2);
  free((void *) u1);
  free((void *) Ah);
  free((void *) p1);
  
    *steps_max=count;

    return(residual);
 }

/*  ==========================================================================  */

 


void v_from_vector(E,V,F)
     struct All_variables *E;
     float **V;
     double *F;
{
  int node,d;
  unsigned int type;

  const int addi_dof = additional_dof[E->mesh.nsd];
  const int nno = E->mesh.nno;
  const int dofs = E->mesh.dof;

  for(node=1;node<=nno;node++)     {
      if(E->node[node] & OFFSIDE) continue;
     
      V[1][node] = F[E->id[node].doff[1]]; 
      V[2][node] = F[E->id[node].doff[2]]; 
      if(dofs==3) V[3][node] = F[E->id[node].doff[3]];
      if (E->node[node] & VBX)
             V[1][node] = E->VB[1][node]; 
      if (E->node[node] & VBZ)
             V[2][node] = E->VB[2][node]; 
      if (dofs==3 && E->node[node] & VBY)
             V[3][node] = E->VB[3][node]; 

    }
  return;
}
