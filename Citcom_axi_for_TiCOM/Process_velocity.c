/*  Here are the routines which process the results of each velocity solution, and call
    the relevant output routines. At this point, the velocity and pressure fields have
    been calculated and stored at the nodes. The only properties of the velocity field
    which are already known are those required to check convergence of the iterative
    scheme and so on. */

#include <math.h>
//#include <malloc.h>
#include <sys/types.h>
//#include <stdlib.h> /* for "system" command */
//#include <malloc.h>
#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#include "element_definitions.h"
#include "global_defs.h"

void process_new_velocity(E,ii)
    struct All_variables *E;
    int ii;
{ 
    void output_velo_related();
    void get_STD_topo();
    void get_CBF_topo();
    void compare_ana();

    static int been_here=0;

    if(been_here==0) {
 	   E->monitor.length_scale = E->data.layer_km/E->mesh.layer[2]; /* km */
	   been_here++;
       } 

    if ( ((ii % E->control.record_every) == 0))     {

      get_STD_topo(E,E->slice.tpg,E->slice.tpgb,ii); 

 /*        get_CBF_topo(E,E->slice.tpg,E->slice.tpgb);   */

      output_velo_related(E,ii);         /* also topo */


/*
      compare_ana(E);
*/

      }

    return;
}

 void compare_ana(E)
  struct All_variables *E;
 {

  int ll1,ll,mm,m,i,it,j,k,snode,node,lev,nox,noy,noz;
  float vt,vf,vt0,vf0,tot,totb,etot,etotb,x,*power[100],z,th,fi,*TG;
  FILE *fp1;
  char output_file[255];

  double y2_surf,y2_botm,y3_surf,y3_botm,aa,con,temp,t1,f1,r1,sphere_h(),Ra1, Ra;
  double modified_plgndr_a(),multis();

  fp1 = fopen ("junk","w");

/* for 2 0, delta funtion source at r=0.75. ri=0.5 and ro=1. constant visc */
  y2_surf = -1.14111e-02;
  y2_botm = 1.36282e-02;
  y3_surf = -0.407752;
  y3_botm = 0.824252;

        mm = E->convection.perturb_mm[0];
        ll = E->convection.perturb_ll[0];

    ll1 = ll + 1;
       aa = 0.0;
       if (mm==0) aa = 1.0;
       Ra=sqrt((2.0-aa)*(2*ll+1)*multis(ll-mm)/(4.0*M_PI*multis(ll+mm)));
       Ra1=sqrt((2.0-aa)*(2*ll1+1)*multis(ll1-mm)/(4.0*M_PI*multis(ll1+mm)));

       for(node=1;node<=E->mesh.nno;node++)   {
          t1=E->X[1][node];
          vt= 0.0;
          if (node%E->mesh.noz==0)  {
             vt = y2_surf*( Ra/Ra1*(ll-mm+1)*modified_plgndr_a(ll+1,mm,t1)
                     - (ll+1)*cos(t1)*modified_plgndr_a(ll,mm,t1) )/(sin(t1)+1e-32);
             }
          else if((node-1)%E->mesh.noz==0)  {
             vt = y2_botm*( Ra/Ra1*(ll-mm+1)*modified_plgndr_a(ll+1,mm,t1)
                     - (ll+1)*cos(t1)*modified_plgndr_a(ll,mm,t1) )/(sin(t1)+1e-32);
             }
         i = node;
        fprintf(fp1,"%5d %.5e %.5e %.6e %.6e %.6e\n",node,E->X[1][i],E->X[2][i],E->V[1][i],vt,E->T[i]); 
          }
       for(i=1;i<=E->mesh.nsf;i++)   {
          node = E->surf_node[i];
          t1 = E->X[1][node];
          vt = -y3_surf*modified_plgndr_a(ll,mm,t1); 
          vf = y3_botm*modified_plgndr_a(ll,mm,t1); 
          if (i==1)  {
              vt0 = vt;
              vf0 = vf;
              }
          fprintf(fp1,"%4d %.4e %.4e %.4e %.4e %.4e\n",i,E->X[1][node],vt-vt0,vf-vf0,E->slice.tpg[i]-E->slice.tpg[1],E->slice.tpgb[i]-E->slice.tpgb[1]);
          }

    fclose(fp1);


    return;
 }
 
/* ===============================================   */

void get_surface_velo(E, SV)
  struct All_variables *E;
  float *SV;
  {

  int el,els,i,m,node,lev;
  char output_file[255];
  FILE *fp;

  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];
  const int nno=E->mesh.nno;

  lev = E->mesh.levmax;

  m = 0;

  for (node=1;node<=nno;node++)
    if ((node-1)%E->mesh.noz==0)   {
      i = (node-1)/E->mesh.noz + 1;
        SV[(i-1)*2+1] = E->V[1][node];
        SV[(i-1)*2+2] = E->V[3][node];
      }

  return;
  }

/* ===============================================   */

void get_ele_visc(E, EV)
  struct All_variables *E;
  float *EV;
  {

  int el,j,lev;

  const int nel=E->mesh.nel;
  const int vpts=vpoints[E->mesh.nsd];

  lev = E->mesh.levmax;

  for (el=1;el<=nel;el++)   {
    EV[el] = 0.0;
    for (j=1;j<=vpts;j++)
      EV[el] +=  E->EVI[lev][(el-1)*vpts+j];

    EV[el] /= vpts;
    }

  return;
  }


void get_surf_stress(E,SXX,SYY,SZZ,SXY,SXZ,SZY)
  struct All_variables *E;
  float *SXX,*SYY,*SZZ,*SXY,*SXZ,*SZY;
  {
  int i,node,stride;

  stride = E->mesh.nsf*6;

  for (node=1;node<=E->mesh.nno;node++)
     if ( ((node-1)%E->mesh.noz)==0 )  {
        i = (node-1)/E->mesh.noz+1;
        E->stress[(i-1)*6+1] = SXX[node];
        E->stress[(i-1)*6+2] = SZZ[node];
        E->stress[(i-1)*6+3] = SYY[node];
        E->stress[(i-1)*6+4] = SXY[node];
        E->stress[(i-1)*6+5] = SXZ[node];
        E->stress[(i-1)*6+6] = SZY[node];
        }
     else if ( ((node-2)%E->mesh.noz)==0 )  {
        i = (node-2)/E->mesh.noz+1;
        E->stress[stride+(i-1)*6+1] = SXX[node];
        E->stress[stride+(i-1)*6+2] = SZZ[node];
        E->stress[stride+(i-1)*6+3] = SYY[node];
        E->stress[stride+(i-1)*6+4] = SXY[node];
        E->stress[stride+(i-1)*6+5] = SXZ[node];
        E->stress[stride+(i-1)*6+6] = SZY[node];
        }

  return;
  }
