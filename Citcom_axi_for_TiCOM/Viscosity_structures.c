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
/* Functions relating to the determination of viscosity field either
   as a function of the run, as an initial condition or as specified from
   a previous file */


#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif
//#include <malloc.h>
#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#ifndef STRING_H
#define STRING_H
#include <string.h>
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

void get_viscosity_option(E)
     struct All_variables *E;
{
    void viscosity_for_system();
    int input_string();
    int input_boolean();
    int input_int();
 
    /* general, essential default */
  
    E->viscosity.update_allowed = 0; 
    E->viscosity.SDEPV = E->viscosity.TDEPV = E->viscosity.CHEMDEPV = 0;
    E->viscosity.EXPX=0;
  
    input_string("Viscosity",E->viscosity.STRUCTURE,NULL);   /* Which form of viscosity */
    
    input_boolean("VISC_EQUIVDD",&(E->viscosity.EQUIVDD),"off");    /* Whether to average it */
    input_int("equivdd_opt",&(E->viscosity.equivddopt),"1");
    input_int("equivdd_x",&(E->viscosity.proflocx),"1");
    input_int("equivdd_y",&(E->viscosity.proflocy),"1");
  
    input_boolean("VISC_SMOOTH",&(E->viscosity.SMOOTH),"off");
    input_int ("visc_smooth_cycles",&(E->viscosity.smooth_cycles),"0");
    
    if ( strcmp(E->viscosity.STRUCTURE,"system") == 0) /* Interpret */ {
      fprintf(E->fp,"Viscosity derived from system state\n");
      E->viscosity.FROM_SYSTEM = 1;
      viscosity_for_system(E);
    }

  return;

}

/* ============================================ */


void viscosity_for_system(E)
     struct All_variables *E;
{
    void get_system_viscosity();
    void twiddle_thumbs();
    int input_int();
    int input_float_vector();
    int input_float();
    int input_float_vector();
    int input_boolean();
    int input_string();
    int l,i;
  
  /* default values .... */

   for(i=0;i<40;i++) {
       E->viscosity.N0[i]=1.0;
       E->viscosity.T[i] = 0.0;
       E->viscosity.Z[i] = 0.0;
       E->viscosity.E[i] = 0.0;
       E->viscosity.T0[i] = 0.0;
   }

  /* read in information */
    input_int("rheol",&(E->viscosity.RHEOL),"essential");
    input_int("num_mat",&(E->viscosity.num_mat),"1");
 
    input_float_vector("viscT",E->viscosity.num_mat,(E->viscosity.T));  /* redundant */
    input_float_vector("viscT1",E->viscosity.num_mat,(E->viscosity.T));
    input_float_vector("viscZ",E->viscosity.num_mat,(E->viscosity.Z));
    input_float_vector("viscE",E->viscosity.num_mat,(E->viscosity.E));
    input_float_vector("viscT0",E->viscosity.num_mat,(E->viscosity.T0));
    input_float_vector("visc0",E->viscosity.num_mat,(E->viscosity.N0)); /* redundant */
    input_float_vector("viscN0",E->viscosity.num_mat,(E->viscosity.N0));
  
    input_boolean("CHEMDEPV",&(E->viscosity.CHEMDEPV),"on");
    input_boolean("TDEPV",&(E->viscosity.TDEPV),"on");
    input_boolean("SDEPV",&(E->viscosity.SDEPV),"off");
    
    input_float("sdepv_misfit",&(E->viscosity.sdepv_misfit),"0.001");
    input_float_vector("sdepv_expt",E->viscosity.num_mat,(E->viscosity.sdepv_expt));
    input_float_vector("sdepv_trns",E->viscosity.num_mat,(E->viscosity.sdepv_trns));
     
    input_boolean("TDEPV_AVE",&(E->viscosity.TDEPV_AVE),"off");
    input_boolean("VFREEZE",&(E->viscosity.FREEZE),"off");
    input_boolean("VMAX",&(E->viscosity.MAX),"off");
    input_boolean("VMIN",&(E->viscosity.MIN),"off");
    input_boolean("VISC_UPDATE",&(E->viscosity.update_allowed),"on");

    if (E->viscosity.TDEPV || E->viscosity.SDEPV)
        E->viscosity.update_allowed = 1;
  
    input_float("freeze_thresh",&(E->viscosity.freeze_thresh),"0.0");
    input_float("freeze_value",&(E->viscosity.freeze_value),"1.0");
    input_float("visc_max",&(E->viscosity.max_value),"nodefault");
    input_float("visc_min",&(E->viscosity.min_value),"nodefault");

    input_boolean("VISC_GUESS",&(E->viscosity.guess),"off");
    input_string("visc_old_file",E->viscosity.old_file," ");

 if (E->rad_heat.num)   {
    for (l=0;l<E->viscosity.num_mat;l++)  {
      E->viscosity.Z[l] = E->data.density*E->data.grav_acc*E->sphere.ro_dim*E->viscosity.Z[l]/(E->data.gas_const*E->data.ref_temperature);
      E->viscosity.E[l] = E->viscosity.E[l]/(E->data.gas_const*E->data.ref_temperature);
     fprintf(E->fp,"E & Z %g %g\n",E->viscosity.E[l],E->viscosity.Z[l]);
     }
   }

    if(!E->viscosity.update_allowed)  {
      get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
      }
 
return;
}


void get_system_viscosity(E,propogate,evisc,visc)
     struct All_variables *E;
     int propogate;
     float *visc,*evisc;     
{ 
    void visc_from_mat();
    void visc_from_T();
    void visc_from_S();
    void apply_viscosity_smoother();
    void v_to_nodes();

    int i,j;
    float *visc_old,*evisc_old;

    const int vpts = vpoints[E->mesh.nsd];

    if(E->viscosity.TDEPV)
       visc_from_T(E,visc,evisc,propogate); 
    else
       visc_from_mat(E,visc,evisc);   

    if(E->viscosity.SDEPV)
       visc_from_S(E,visc,evisc,propogate);
    
/*    if (E->viscosity.SMOOTH) 
       apply_viscosity_smoother(E,visc,evisc);
*/
    if(E->viscosity.MAX) {
      for(i=1;i<=E->mesh.nel;i++)
          for(j=1;j<=vpts;j++) {
            if(evisc[(i-1)*vpts + j] > E->viscosity.max_value)
               evisc[(i-1)*vpts + j] = E->viscosity.max_value;
            }
      }

    if(E->viscosity.MIN) {
      for(i=1;i<=E->mesh.nel;i++)
          for(j=1;j<=vpts;j++)
            if(evisc[(i-1)*vpts + j] < E->viscosity.min_value)
               evisc[(i-1)*vpts + j] = E->viscosity.min_value;
      }

    v_to_nodes(E,evisc,visc,E->mesh.levmax);

 return;
}


void apply_viscosity_smoother(E,visc,evisc)
     struct All_variables *E;
     float *visc,*evisc;

{
    void p_to_centres();
    void p_to_nodes();
    
    float  *ViscCentre;
    int i;

    ViscCentre = (float *)malloc((E->mesh.nno+10)*sizeof(float));
   
  
    for(i=1;i<=E->viscosity.smooth_cycles;i++)  {
	p_to_centres(E,visc,ViscCentre,E->mesh.levmax);
	p_to_nodes(E,ViscCentre,visc,E->mesh.levmax);
    }

 
    free ((void *)ViscCentre);
 
  return;
}

void visc_from_mat(E,Eta,EEta)
     struct All_variables *E;
     float *Eta,*EEta;
{

    int i,j,k,l,z,jj,kk;
   float temp1,temp[9],C1,C2,visc1,pres0,slope,temp2,rii[9];
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    for(i=1;i<=E->mesh.nel;i++)   {
      for(jj=1;jj<=vpoints[E->mesh.nsd];jj++)  {
        EEta[ (i-1)*vpoints[E->mesh.nsd]+jj ]=E->viscosity.N0[E->mat[i]-1];
        }
      }

    return;
}

void visc_from_T(E,Eta,EEta,propogate)
     struct All_variables *E;
     float *Eta,*EEta;
     int propogate;

{
    void remove_horiz_ave();
    int i,j,k,e,l,z,jj,kk,imark;
    float c1,c2,c3,zero,e_6,one,eta0,slope,Tave,depth,temp[9],rii[9],tempa,TT[9];
    double temp1,temp2,temp3;
    static int visits=0;
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->mesh.nel;
    const int noz = E->mesh.noz;

    one = 1.0;
    zero = 0.0;

    slope = (1.0-E->data.visc_factor)/(E->sphere.ro-E->sphere.ri);

    for(i=noz;i>=1;i--)
       E->Have.Tadi[i] = 0.0;

    E->data.T_adi0 = 0;
    E->data.T_adi1 = 0;

    switch (E->viscosity.RHEOL)   {
    case 1:
        if(propogate && visits==0) {
            fprintf(E->fp,"\tRheological option 1:\n");

            for(l=1;l<=E->viscosity.num_mat;l++) {
              fprintf(E->fp,"\tlayer %d/%d: E=%g T1=%g \n",
                      l,E->viscosity.num_mat,
                      E->viscosity.E[l-1],E->viscosity.T[l-1]);
            }
            fflush(E->fp);
        }

     temp2 = one;
     if (E->control.adi_heating)  {
	  temp3 = 0;
	  E->data.T_adi0 = 0;
	  for(i=noz;i>1;i--)
            if (E->X[2][i]<(2*E->viscosity.zlith-E->sphere.ro))  {
               temp3 = E->Have.T[i];
               break;
               }

          E->data.T_adi0 = temp3;
	  E->Have.Tadi[noz] = temp3;
	  for(i=noz;i>1;i--)   {
            if (E->X[2][i]<(2*E->viscosity.zlith-E->sphere.ro))
               temp3 = temp3 + E->data.disptn_number*(E->expansivity[i]+E->expansivity[i-1])*(E->Have.T[i]+E->Have.T[i-1]+2*E->data.surf_temp)/4.0*(E->X[2][i]-E->X[2][i-1]);
            E->Have.Tadi[i-1] = temp3;
	    }

	  E->data.T_adi1 = E->Have.Tadi[1];
	  temp2 = one-E->data.T_adi1 + E->data.T_adi0;
     }


        for(i=1;i<=nel;i++)   {
            l = E->mat[i];
	    e = (i-1)%E->mesh.elz+1;

            tempa = E->viscosity.N0[l-1];
	    temp1 = (E->Have.Tadi[e] + E->Have.Tadi[e+1])*0.5 -E->data.T_adi0;
	    temp1 = 0;

            for(kk=1;kk<=ends;kk++)
               TT[kk] = E->T[E->ien[i].node[kk]];

	    for(jj=1;jj<=vpts;jj++) {
		       temp[jj] = 0;
		       rii[jj] = 0;
		       for(kk=1;kk<=ends;kk++)  {
			  temp[jj] += E->T[E->ien[i].node[kk]]*E->N.vpt[GNVINDEX(kk,jj)];
                          rii[jj] += E->X[2][E->ien[i].node[kk]]*E->N.vpt[GNVINDEX(kk,jj)];
			  }
		       EEta[ (i-1)*vpoints[E->mesh.nsd]+jj ]=tempa*exp(E->viscosity.E[l-1]*(temp2-(temp[jj]))/temp2)
			       *(E->data.visc_factor+slope*(rii[jj]-E->sphere.ri));
		       }
		   }
		break;

    case 2:
        fprintf(stderr,"not supported for rheol=2\n");
    	break;

    case 3:
        if(propogate && visits==0) {
            fprintf(E->fp,"\tRheological option 3:\n");

            for(l=1;l<=E->viscosity.num_mat;l++) {
              fprintf(E->fp,"\tlayer %d/%d: E=%g T1=%g \n",
                      l,E->viscosity.num_mat,
                      E->viscosity.E[l-1],E->viscosity.T[l-1]);
            }
            fflush(E->fp);
        }

	/*fprintf(stderr,"Ts=%g\n",E->data.surf_temp);*/
	
        for(i=1;i<=nel;i++)   {
            l = E->mat[i];
            tempa = E->viscosity.N0[l-1];
            j = 0;

            for(kk=1;kk<=ends;kk++)
               TT[kk] = E->T[E->ien[i].node[kk]];

            for(jj=1;jj<=vpts;jj++) {
               rii[jj] = 0;
               temp[jj] = 0;
               for(kk=1;kk<=ends;kk++)  {
                  rii[jj] += E->X[2][E->ien[i].node[kk]]*E->N.vpt[GNVINDEX(kk,jj)];
                  temp[jj] += E->T[E->ien[i].node[kk]]*E->N.vpt[GNVINDEX(kk,jj)];
                  }
               temp1 = (E->viscosity.E[E->mat[i]-1]+ 
                   E->viscosity.Z[E->mat[i]-1]*(E->sphere.ro-rii[jj]))
                          /(temp[jj]+E->data.surf_temp)
                 -(E->viscosity.E[E->mat[i]-1]+
                   E->viscosity.Z[E->mat[i]-1]*(E->sphere.ro-E->sphere.ri))
                          /(1.0+E->data.surf_temp);
               EEta[ (i-1)*vpoints[E->mesh.nsd]+jj ]=E->viscosity.N0[E->mat[i]-1]*exp(temp1);
               }
           
           }
        break;
	
    }

    visits++;

  return;  }


void visc_from_S(E,Eta,EEta,propogate)
     struct All_variables *E;
     float *Eta,*EEta;
     int propogate;
{
    static int visits = 0;
    float one,two,scale,stress_magnitude,depth,exponent1;
    float *eedot;
    
    void strain_rate_2_inv();
    int e,l,z,jj,kk;

    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->mesh.nel;

    eedot = (float *) malloc((2+nel)*sizeof(float));
    one = 1.0;
    two = 2.0;
    
	  if (visits==0)   {
        for(e=1;e<=nel;e++)
            eedot[e] = one; 
        } 
      else
        strain_rate_2_inv(E,eedot,1);

      for(e=1;e<=nel;e++)   {
        exponent1= one-one/E->viscosity.sdepv_expt[E->mat[e]-1];
	    scale=pow(two*eedot[e]/E->viscosity.sdepv_trns[E->mat[e]-1],exponent1);
        for(jj=1;jj<=vpts;jj++)
	      EEta[(e-1)*vpts + jj] = two*EEta[(e-1)*vpts+jj]/(one+scale*pow(EEta[(e-1)*vpts+jj],exponent1) );
	    }


       visits ++;

	free ((void *)eedot);
    return;  
}



void strain_rate_2_inv(E,EEDOT,SQRT)
     struct All_variables *E;
     float *EEDOT;
     int SQRT;
{
    void get_global_shape_fn();
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    
    double aaa,xk[3][5],edot[4][4],dudx[4][4];
    float VV[4][9];

    int e,i,p,q,n,nel,k;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->mesh.nno;
    const int vpts = vpoints[dims];
 
    nel = E->mesh.nel;

    for(e=1;e<=nel;e++) {
  
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,xk,2,lev);

      for(i=1;i<=ends;i++)   {
        n=E->ien[e].node[i];
          VV[1][i] = E->V[1][n];
          VV[2][i] = E->V[2][n];
        }
      
      dudx[1][1] = 0.0;  
      dudx[1][2] = 0.0;  
      dudx[2][2] = 0.0;  
      dudx[3][3] = 0.0;  
      
      for(i=1;i<=ends;i++)  {
        dudx[1][1] += VV[1][i]*GNx.ppt[GNPXINDEX(0,i,1)] + VV[2][i]*E->N.ppt[GNPINDEX(i,1)]*xk[2][1];
        dudx[2][2] += VV[2][i]*GNx.ppt[GNPXINDEX(1,i,1)];
        dudx[3][3] += ( VV[2][i] + VV[1][i]*xk[1][1] )*E->N.ppt[GNPINDEX(i,1)]*xk[2][1];
        dudx[1][2] += VV[2][i]*GNx.ppt[GNPXINDEX(0,i,1)] + VV[1][i]*GNx.ppt[GNPXINDEX(1,i,1)]-VV[1][i]*xk[2][1]*E->N.ppt[GNPINDEX(i,1)];
        }

      edot[1][1] = 2.0*dudx[1][1];  
      edot[2][2] = 2.0*dudx[2][2];  
      edot[3][3] = 2.0*dudx[3][3];  
      edot[1][2] = dudx[1][2];  
	

         EEDOT[e] = edot[1][1]*edot[1][1] + edot[2][2]*edot[2][2]
                  + edot[3][3]*edot[3][3] + edot[1][2]*edot[1][2]*2.0; 

      }

    if(SQRT)
	for(e=1;e<=nel;e++)
	    EEDOT[e] =  sqrt(0.5 *EEDOT[e]);
    else
	for(e=1;e<=nel;e++)
	    EEDOT[e] *=  0.5;
  
    return;
}

 int layers(E,x2)
    struct All_variables *E;
    float x2;
 {
    int llayers = 0;

    if (x2>=E->viscosity.zlith)
      llayers = 3;
    else if (x2>=E->viscosity.zlm)
      llayers = 2;

    return (llayers);
    }

