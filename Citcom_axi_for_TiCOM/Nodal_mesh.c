/* Functions relating to the building and use of mesh locations ... */


#include <math.h>
//#include <malloc.h>
#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

/* =================================================
   Standard node positions including mesh refinement 

   =================================================  */

void node_locations(E)
     struct All_variables *E;
{ 
  int lev,i,j,k,ijk[4],ii,d,node;
  float *XX[4],*XG[4],dx[4],dxx[40],dx1,dx2,dc,rc,dr1,dr2,dr3,dr4,dr0,dr5,dr6;
  int nox,noz,noy,fn,step,ncr;

  const int dims = E->mesh.nsd;

  int input_int();
  int input_float_vector();
  int input_int_vector();

  input_int("r_grid_layers",&(E->segment.zlayers),"1");
  input_float_vector("rr",E->segment.zlayers,(E->segment.zzlayer));
  input_int_vector("nr",E->segment.zlayers,(E->segment.nzlayer));

  input_int("theta_grid_layers",&(E->segment.xlayers),"1");
  input_float_vector("theta",E->segment.xlayers,(E->segment.xxlayer));
  input_int_vector("ntheta",E->segment.xlayers,(E->segment.nxlayer));


     for(d=1;d<=E->mesh.nsd;d++) {
       XX[d] = (float *)malloc((2+E->mesh.nnx[d])*sizeof(float)); 
       XG[d] = (float *)malloc ((2+1)*sizeof(float)); 
       }

     dx[1] = E->mesh.layer[1]/(E->mesh.nnx[1]-1);
     XX[1][1] = 0.0;
     for(i=2;i<=E->mesh.nnx[1];i++)
  	      XX[1][i] = XX[1][i-1]+dx[1];

     dx[2] = (E->sphere.ro-E->sphere.ri)/(E->mesh.nnx[2]-1);
     XX[2][1] = E->sphere.ri;
     for(i=2;i<=E->mesh.nnx[2];i++)
  	      XX[2][i] = XX[2][i-1]+dx[2];

  for (j=1;j<E->segment.xlayers;j++)
    dxx[j] = (E->segment.xxlayer[j]-E->segment.xxlayer[j-1])*M_PI
            /(E->segment.nxlayer[j]-E->segment.nxlayer[j-1]);
  j=1;
  for(i=2;i<E->mesh.nnx[1];i++)   {
    if (i<=E->segment.nxlayer[j])
       XX[1][i] = XX[1][i-1] + dxx[j];
    if (i==E->segment.nxlayer[j])
       j++;
    }

  for (j=1;j<E->segment.zlayers;j++)
    dxx[j] = (E->segment.zzlayer[j]-E->segment.zzlayer[j-1])
            /(E->segment.nzlayer[j]-E->segment.nzlayer[j-1]);
  j=1;
  for(i=2;i<E->mesh.nnx[2];i++)   {
    if (i<=E->segment.nzlayer[j])
       XX[2][i] = XX[2][i-1] + dxx[j];
    if (i==E->segment.nzlayer[j])
       j++;
    }

   for(d=1;d<=E->mesh.nsd;d++)
     for(i=1;i<=E->mesh.nnx[d];i++)
       E->XP[d][i] = XX[d][i]; 
    
   for (lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {

     nox=E->mesh.NOX[lev]; 
     noy=E->mesh.NOY[lev];
     noz=E->mesh.NOZ[lev];

    if (E->control.NMULTIGRID||E->control.EMULTIGRID)
        step = (int) pow(2.0,(double)(E->mesh.levmax-lev));
    else
        step = 1;

     for(ijk[1]=1;ijk[1]<=nox;ijk[1]++)
       for(ijk[2]=1;ijk[2]<=noz;ijk[2]++)
         for(ijk[3]=1;ijk[3]<=noy;ijk[3]++)   {
           node=ijk[2]+(ijk[1]-1)*noz+(ijk[3]-1)*noz*nox;
           for(d=1;d<=E->mesh.nsd;d++)
             E->XX[lev][d][node] = XX[d][(ijk[d]-1)*step+1]; 
           }
      }

     for(d=1;d<=E->mesh.nsd;d++) {
       free((void *)XX[d]);
       free((void *)XG[d]);
       }

  if (E->control.verbose) 
    for (lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {
      fprintf(E->fp,"output_coordinates %d\n",lev);
      if (dims==2)  {
         for (i=1;i<=E->mesh.NNO[lev];i++)
             fprintf(E->fp,"%d %g %g\n",i,E->XX[lev][1][i],E->XX[lev][2][i]);
         }
      else if (dims==3)  {
         for (i=1;i<=E->mesh.NNO[lev];i++)
             fprintf(E->fp,"%d %g %g %g\n",i,E->XX[lev][1][i],E->XX[lev][2][i],E->XX[lev][3][i]);
         }
      }

return;  }


/********************************************************
 * update_node_locations                                *
 *      function to determine the new radial position   *
 *      of the nodes after a change in thickness of     *
 *      the model domain. Horizontal positions do not   *
 *      change. The old node locations are NOT          *
 *      overwritten at this point, because we still     *
 *      need to interpolate and we need both sets of    *
 *      coordinates for that.                           *
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *                                                      *
 * Returns                                              *
 *      none                                            *
 ********************************************************/
void update_node_locations(E)
     struct All_variables *E;
{ 
  int lev,e,i,j,k,ijk[4],ii,d,node;
  float *RR,RR1,dr,drr[40];
  int nox,noz,noy,step;

  const int dims = E->mesh.nsd;

  /* Allocate new X array in r-direction only */
  RR = (float *)malloc((2+E->mesh.nnx[2])*sizeof(float)); 

  /* find radial positions of all nodes */
  fprintf(stderr,"ri %f -",E->sphere.ri);
  E->sphere.ri -= E->sphere.deltarb; /* Adjust inner radius */
  fprintf(stderr,"dh %f =",E->sphere.deltarb);
  E->sphere.ri = max(E->sphere.ri,E->sphere.rcore); /* Can't be below core */
  fprintf(stderr," %f rc %f\n",E->sphere.ri,E->sphere.rcore);
  //dr = (E->sphere.ro-E->sphere.ri+E->sphere.deltarb)/(E->mesh.nnx[2]-1);
  dr = (E->sphere.ro-E->sphere.ri)/(E->mesh.nnx[2]-1);
  fprintf(stderr,"dr %f\n",dr);
  //RR[1] = max(E->sphere.rcore,(E->sphere.ri-E->sphere.deltarb));
  RR[1] = E->sphere.ri;
//  fprintf(stderr," %f %f\n",E->sphere.ri,RR[1]);
//  fprintf(stderr,"i %d rold %f rnew %f\n",1,(E->X[2][1]),(RR[1]));
  for(i=2;i<=E->mesh.nnx[2];i++){
      RR[i] = RR[i-1]+dr;
//      fprintf(stderr,"i %d RRi %f RRi-1 %f dr %f\n",i,RR[i],RR[i-1],dr);
//      fprintf(stderr,"i %d rold %f rnew %f\n",i,(E->X[2][i]),(RR[i]));
      }

  /* Adjust bottom grid layer */
  E->segment.zzlayer[0] -= E->sphere.deltarb;
  for (j=1;j<E->segment.zlayers;j++)
    drr[j] = (E->segment.zzlayer[j]-E->segment.zzlayer[j-1])
            /(E->segment.nzlayer[j]-E->segment.nzlayer[j-1]);
  j=1;
  for(i=2;i<E->mesh.nnx[2];i++)   {
    if (i<=E->segment.nzlayer[j]) {
       RR[i] = RR[i-1] + drr[j];
      }
    if (i==E->segment.nzlayer[j])
       j++;
  }

//  for(i=1;i<=E->mesh.nnx[2];i++)
//      fprintf(stderr,"i %d rold %f rnew %f\n",i,(E->X[2][i]),(RR[i]));

  /* Copy old E->X[2] at top mg level into a backup array */
  for(i=1;i<=E->mesh.nno;i++)
    E->newvars.oldX[i] = E->X[2][i]; /* Radial coord only */
  /* Do the same for old E->ECO.centre[2] at top mg level */
 for (e=1;e<=E->mesh.nel;e++)
    E->newvars.oldecocentre[e] = E->eco[e].centre[2]; /* Radial coord only */
   
  /* New XP only in radial direction*/
  for(i=1;i<=E->mesh.nnx[2];i++)
    E->newvars.XP[i] = RR[i]; 

  /* Define new XX for radial coord only at all mg levels */
  for (lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {
  
    if (E->control.NMULTIGRID||E->control.EMULTIGRID)
      step = (int) pow(2.0,(double)(E->mesh.levmax-lev));
    else
      step = 1;

    nox=E->mesh.NOX[lev]; 
    noy=E->mesh.NOY[lev];
    noz=E->mesh.NOZ[lev];
    
    for(ijk[1]=1;ijk[1]<=nox;ijk[1]++)
      for(ijk[2]=1;ijk[2]<=noz;ijk[2]++)
        for(ijk[3]=1;ijk[3]<=noy;ijk[3]++)   {
          node=ijk[2]+(ijk[1]-1)*noz+(ijk[3]-1)*noz*nox;
          
          E->newvars.XX[lev][node] = RR[(ijk[2]-1)*step+1]; 
          E->XX[lev][2][node] = RR[(ijk[2]-1)*step+1]; 
        }
  }

  E->newvars.coords = 1;              /* Set flag to output coords again */

  if (E->control.verbose) 
    for (lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {
      fprintf(E->fp,"output_new_coordinates %d\n",lev);
      if (dims>=1)  {
         for (i=1;i<=E->mesh.NNO[lev];i++)
             //fprintf(E->fp,"%d %g\n",i,E->newvars.XX[lev][i]);
             fprintf(E->fp,"%d %g %g\n",i,E->XX[lev][1][i],E->XX[lev][2][i]);
      }
    }

  free((void *)RR);

  return;
}

void flogical_mesh_to_real(E,data,level)
     struct All_variables *E;
     float *data;
     int level;

{ int i,j,n1,n2;

  return;
}

void dp_to_nodes(E,P,PN,lev)
     struct All_variables *E;
     double *P;
     float *PN;
     int lev;

{ int e,element,node,j;

  for(node=1;node<=E->mesh.NNO[lev];node++)
    PN[node] =  0.0;
	  
  for(element=1;element<=E->mesh.NEL[lev];element++) {

      for(j=1;j<=enodes[E->mesh.nsd];j++)  {
     	  node = E->IEN[lev][element].node[j];
    	  PN[node] += P[element] * E->TW[lev][node] ; 
    	  }

      } 

     return; }

void p_to_nodes(E,P,PN,lev)
     struct All_variables *E;
     float *P,*PN;
     int lev;

{ int e,element,node,j;

  for(node=1;node<=E->mesh.NNO[lev];node++)
    PN[node] =  0.0;
	  
  for(element=1;element<=E->mesh.NEL[lev];element++) {

      for(j=1;j<=enodes[E->mesh.nsd];j++)  {
     	  node = E->IEN[lev][element].node[j];
    	  PN[node] += P[element] * E->TW[lev][node] ; 
    	  }

      } 

     return; }


void p_to_centres(E,PN,P,lev)
     struct All_variables *E;
     float *PN,*P;
     int lev;

{  int p,element,node,j;
   double weight;

   for(p=1;p<=E->mesh.NEL[lev];p++)
     P[p] = 0.0;

   weight=1.0/((double)enodes[E->mesh.nsd]) ;
   
   for(p=1;p<=E->mesh.NEL[lev];p++)
     for(j=1;j<=enodes[E->mesh.nsd];j++)
       P[p] +=  PN[E->IEN[lev][p].node[j]] * weight;

   return;  
   }


void v_to_intpts(E,VN,VE,lev)
  struct All_variables *E;
  float *VN,*VE;
  int lev;
  {

   int e,i,j,k;
   const int nsd=E->mesh.nsd;
   const int vpts=vpoints[nsd];
   const int ends=enodes[nsd];

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)                 {
        VE[(e-1)*vpts + i] = 0.0;
        for(j=1;j<=ends;j++)
          VE[(e-1)*vpts + i] += VN[E->IEN[lev][e].node[j]] *  E->N.vpt[GNVINDEX(j,i)];
        }

   return;
  }

void v_to_nodes(E,VE,VN,lev)
   struct All_variables *E;
   float *VE,*VN;
   int lev;
   {
    int e,i,j,k,n;
    const int nsd=E->mesh.nsd;
    const int vpts=vpoints[nsd];
    const int ends=enodes[nsd];

    for(i=1;i<=E->mesh.NNO[lev];i++)
	    VN[i] = 0.0;

    for(e=1;e<=E->mesh.NEL[lev];e++)
      for(j=1;j<=ends;j++) {
        n = E->IEN[lev][e].node[j];
        for(i=1;i<=vpts;i++)
          VN[n] += E->N.vpt[GNVINDEX(j,i)] * E->TW[lev][n] * VE[(e-1)*vpts + i];
        }

    return;
    }

void visc_to_intpts(E,VN,VE,lev)
   struct All_variables *E;
   float *VN,*VE;
   int lev;
   {

   int e,i,j,k;
   const int nsd=E->mesh.nsd;
   const int vpts=vpoints[nsd];
   const int ends=enodes[nsd];

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++) {
        VE[(e-1)*vpts + i] = 0.0;
	for(j=1;j<=ends;j++)
          VE[(e-1)*vpts + i] += log(VN[E->IEN[lev][e].node[j]]) *  E->N.vpt[GNVINDEX(j,i)];
        VE[(e-1)*vpts + i] = exp(VE[(e-1)*vpts + i]);
        }

  }


void visc_to_nodes(E,VE,VN,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {
  int e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

  for(i=1;i<=E->mesh.NNO[lev];i++)
    VN[i] = 0.0;

  for(e=1;e<=E->mesh.NEL[lev];e++)
    for(j=1;j<=ends;j++) {
      n = E->IEN[lev][e].node[j];
      temp_visc=0.0;
      for(i=1;i<=vpts;i++)
	temp_visc += E->TW[lev][n] * log(E->N.vpt[GNVINDEX(j,i)] * VE[(e-1)*vpts + i]);
      VN[n] += exp(temp_visc);
      }
   return;
}

void visc_from_gint_to_nodes(E,VE,VN,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {
  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(i=1;i<=E->mesh.NNO[lev];i++)
     VN[i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(j=1;j<=ends;j++)                {
       n = E->IEN[lev][e].node[j];
       temp_visc=0.0;
       for(i=1;i<=vpts;i++)
         temp_visc += E->TW[lev][n] * E->N.vpt[GNVINDEX(j,i)] * VE[(e-1)*vpts + i];
       VN[n] += temp_visc;
       }

   return;
}


 void visc_from_nodes_to_gint(E,VN,VE,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {

  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;
   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)
       VE[(e-1)*vpts+i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)      {
       temp_visc=0.0;
       for(j=1;j<=ends;j++)
         temp_visc += E->N.vpt[GNVINDEX(j,i)]*VN[E->IEN[lev][e].node[j]];

       VE[(e-1)*vpts+i] = temp_visc;
       }

   return;
   }

void visc_from_gint_to_ele(E,VE,VN,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {
  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(i=1;i<=E->mesh.NEL[lev];i++)
     VN[i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)   {
     temp_visc=0.0;
     for(i=1;i<=vpts;i++)
        temp_visc += VE[(e-1)*vpts + i];
     temp_visc = temp_visc/vpts;

     VN[e] = temp_visc;
    }

   return;
}


 void visc_from_ele_to_gint(E,VN,VE,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {

  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)
       VE[(e-1)*vpts+i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)      {

       VE[(e-1)*vpts+i] = VN[e];
       }

   return;
 }

