/* Routine to process the output of the finite element cycles 
   and to turn them into a coherent suite  files  */


#include <fcntl.h>
#include <math.h>
//#include <malloc.h>
//#include <malloc.h>
#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif
//#include <stdlib.h>             /* for "system" command */
#ifndef __sunos__               /* string manipulations */
#include <strings.h>
#else
#include <string.h>
#endif

#ifndef __ELEMENT_DEFINITIONS_H__
#define __ELEMENT_DEFINITIONS_H__
#include "element_definitions.h"
#endif

#ifndef __GLOBAL_DEFS_H__
#define __GLOBAL_DEFS_H__
#include "global_defs.h"
#endif

void output_velo_related(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  FILE *fp1, *fp0;
   float surf,botm;
  int el,els,i,j,k,ii,m,node,fd;
  int nox,noz,noy,nfx,nfz,nfy1,nfy2,size1,size2;
  char output_file[255];
  static float *SV,*EV,*logVi;
  static int been_here=0;

  void get_surface_velo ();
  void get_ele_visc ();
  void coordinates_dx(); 
  void frames_for_dx(); 
  void return_horiz_ave();
  void secular_cooling();


  const int nno = E->mesh.nno;

  if (been_here==0) {
    SV = (float *)malloc((nno+1)*sizeof(float));
    logVi = (float *)malloc((nno+1)*sizeof(float));
  }

  if (been_here==0 && E->control.restart==0 || E->newvars.coords) {
    been_here++;
//    E->newvars.coords = 0;

    sprintf(output_file,"%s/coord.%d",E->control.data_file,file_number);
    fp0=fopen(output_file,"w");
    fprintf(fp0,"%6d %6d %.5e\n",E->mesh.nno,E->advection.timesteps,E->monitor.elapsed_time);
    for (i=1;i<=E->mesh.nno;i++)
      fprintf(fp0,"%6d %.3e %.3e %.5e %.5e %.5e %.5e %.4e\n",i,E->X[1][i],E->X[2][i],E->V[1][i],E->V[2][i],E->T[i],E->C[i],E->Vi[i]);
    fclose(fp0);

    sprintf(output_file,"%s/coord_elem.%d",E->control.data_file,file_number);
    fp0=fopen(output_file,"w");
    fprintf(fp0,"%6d %6d %.5e\n",E->mesh.nel,E->advection.timesteps,E->monitor.elapsed_time);
    for (i=1;i<=E->mesh.nel;i++)
      fprintf(fp0,"%6d %.3e %.3e\n",i,E->eco[i].centre[1],E->eco[i].centre[2]);
    fclose(fp0);
    
    }


  if ((E->advection.timesteps%(E->control.record_every)) == 0)   {


    sprintf(output_file,"%s/temp_comp.%d",E->control.data_file,file_number);
    fp0=fopen(output_file,"w");

    fprintf(fp0,"%6d %6d %.5e %.5e %.5e\n",E->mesh.nno,E->advection.timesteps,E->monitor.elapsed_time,E->monitor.deltah);

    if (E->control.composition) {
      if ((E->advection.timesteps%(20*E->control.record_every)) == 0)  {
        for (i=1;i<=E->mesh.nno;i++) {
					k = ((i-1) % E->mesh.noz) + 1;
          fprintf(fp0,"%.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
										E->T[i],E->C[i],E->Fm[i],E->C_prev[i],
										(E->T[i]-E->solidus[k]),E->V[1][i],E->V[2][i]);
					}
        for (i=1;i<=E->mesh.nel;i++)
          fprintf(fp0,"%.4e\n",E->P[i]);
        }
      else
        for (i=1;i<=E->mesh.nno;i++) {
					k = ((i-1) % E->mesh.noz) + 1;
          fprintf(fp0,"%.4e %.4e %.4e %.4e %.4e\n",E->T[i],E->C[i],E->Fm[i],
										E->C_prev[i],(E->T[i]-E->solidus[k]));
				}
      }
    else  {
      if ((E->advection.timesteps%(20*E->control.record_every)) == 0)  {
        for (i=1;i<=E->mesh.nno;i++)
          fprintf(fp0,"%.5e %.4e %.4e %.4e\n",E->T[i],E->heatflux[i],E->V[1][i],E->V[2][i]);
        for (i=1;i<=E->mesh.nel;i++)
          fprintf(fp0,"%.4e\n",E->P[i]);
        }
      else
        for (i=1;i<=E->mesh.nno;i++)
          fprintf(fp0,"%.5e %.4e\n",E->T[i],E->heatflux[i]);
      }

    fclose(fp0);

    }

 if (E->control.composition && strcmp(E->control.comp_adv_method,"particle")==0 && (E->advection.timesteps%(1*E->control.record_every) == 0) && (E->advection.timesteps >= 0) ){

    sprintf(output_file,"%s/traces.%d",E->control.data_file,file_number);
    fp0=fopen(output_file,"w");
    fprintf(fp0,"%6d %6d %.5e\n",E->advection.markers,
						E->advection.timesteps,E->monitor.elapsed_time);
    for (i=1;i<=E->advection.markers;i++)
      /*fprintf(fp0,"%.5e %.5e %d %d\n",E->XMC[1][i],E->XMC[2][i],
							E->CElement[i],E->C12[i]);*/
      fprintf(fp0,"%.5e %.5e %d %g\n",E->XMC[1][i],E->XMC[2][i],
							E->CElement[i],E->C12f[i]);
		/*if(E->advection.timesteps > 0)
    for (i=1;i<=E->mesh.nel;i++)
      fprintf(fp0,"%d %g %d %d\n",i,E->CE[i],E->advection.element[0][i],
							E->advection.element[1][i]);*/
    fclose(fp0);
		}

if (E->control.composition && (strcmp(E->control.comp_adv_method,"particle")==0) && (E->advection.timesteps%(E->control.record_every)) == 0)  {
		if(E->advection.timesteps >= (0)) {
      sprintf(output_file,"%s/etraces.%d",E->control.data_file,
							file_number);
    	fp0=fopen(output_file,"w");
    	fprintf(fp0,"%d %.5e %d %d %d | %d %d\n",
							E->advection.timesteps,E->monitor.elapsed_time,
							E->advection.markers,E->advection.marker_type[0],
							E->advection.marker_type[1],E->advection.marker_type_prev[0],
							E->advection.marker_type_prev[1]);
      for (i=1;i<=E->mesh.nel;i++)
        fprintf(fp0,"%d %g %g %g %g %d %d | %d %d\n",i,E->eco[i].area,
					E->CE[i],E->FmE[i],E->CE_prev[i],E->advection.element[0][i],
					E->advection.element[1][i],E->advection.element_prev[0][i],
					E->advection.element_prev[1][i]);

    	fclose(fp0);
  	}
	}


  if ((E->advection.timesteps%(E->control.record_every)) == 0)   {

    surf=0.0;
    botm=0.0;
    for (i=1;i<=E->mesh.nox;i++)  {
      j=(int) i*E->mesh.noz;
      if (i>1) {
        surf += ((float) E->slice.shflux[i]+(float) E->slice.shflux[i-1])*0.5*
           sin(0.5*((float) E->X[1][j]+(float) E->X[1][j-(int) E->mesh.noz]))*
           (E->X[1][j]-E->X[1][j-(int) E->mesh.noz]);
        botm += (E->slice.bhflux[i]+E->slice.bhflux[i-1])*0.5*
           sin(0.5*(E->X[1][j]+E->X[1][j-E->mesh.noz]))*
           (E->X[1][j]-E->X[1][j-E->mesh.noz]);
        }
      }
    surf = surf/(1.0-cos(E->X[1][E->mesh.nno]));
    botm = botm/(1.0-cos(E->X[1][E->mesh.nno]));

    /* Adjust CMB Temperature based on CMB heat flux */

    if(E->control.secular)
      secular_cooling(E, botm);
    

    for (i=1;i<=E->mesh.nno;i++){
      SV[i] = sqrt(E->V[1][i]*E->V[1][i] + E->V[2][i]*E->V[2][i]);
      logVi[i] = log(E->Vi[i]);
    }

    return_horiz_ave(E,SV,E->Have.vrms);
    return_horiz_ave(E,logVi,E->Have.Vi);
    for (i=1;i<=E->mesh.noz;i++)
			E->Have.Vi[i] = exp(E->Have.Vi[i]);
    return_horiz_ave(E,E->T,E->Have.T);
    return_horiz_ave(E,E->C,E->Have.Rho);
/*    E->Total.melt_prod = 0.0;
    for (i=1;i<=E->mesh.elx;i++)
			E->Total.melt_prod += E->slice.melt[i];
*/
    sprintf(output_file,"%s/ave.%d",E->control.data_file,file_number);
    fp1=fopen(output_file,"w");

    fprintf(fp1,"%6d %6d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
						E->mesh.nno,E->advection.timesteps,E->monitor.elapsed_time,
						surf,botm,E->rad_heat.total,E->data.T_adi0,E->data.T_adi1,
						E->Total.melt_prod,E->Total.bulk_comp,E->Total.bulk_comp_prev,
            E->sphere.deltarb);
    for (i=1;i<=E->mesh.noz;i++)
      fprintf(fp1,"%.4e %.5e %.5e %.5e %.5e %.5e %.4e %.4e\n",E->X[2][i],E->Have.T[i],E->Have.vrms[i],E->Have.f[i],E->Have.F[i],E->Have.Rho[i],E->Have.Vi[i],E->Have.Tadi[i]);
    fclose(fp1);

    sprintf(output_file,"%s/surf.%d",E->control.data_file,file_number);
    fp1=fopen(output_file,"w");
    for (i=1;i<=E->mesh.nox;i++)  {
      j=i*E->mesh.noz;
      fprintf(fp1,"%.4e %.5e %.5e %.5e %.5e\n",E->X[1][j],E->slice.tpg[i],E->slice.tpgb[i],E->slice.shflux[i],E->slice.bhflux[i]);
      }
      
    fclose(fp1);
    sprintf(output_file,"%s/esurf.%d",E->control.data_file,file_number);
    fp1=fopen(output_file,"w");
    for (i=1;i<=E->mesh.elx;i++)  {
      j=i*E->mesh.elz;
      fprintf(fp1,"%.4e %.5e %.5e\n",E->eco[j].centre[1],E->slice.melt[i],E->slice.new_melt[i]);
      }   

    fclose(fp1);
   }


return;
}

/*      ----------------------------------- */
 void coordinates_dx(E)
  struct All_variables *E;
 {

 int i;

  E->ibm_dx.x1 = (float *) malloc((E->mesh.nno+1)*sizeof(float));
  E->ibm_dx.x2 = (float *) malloc((E->mesh.nno+1)*sizeof(float));

  E->ibm_dx.nox = E->mesh.nox;
  E->ibm_dx.noz = E->mesh.noz;

   for (i=1;i<=E->mesh.nno;i++)   {
      E->ibm_dx.x1[i] = E->X[2][i] * sin(E->X[1][i]);
      E->ibm_dx.x2[i] = E->X[2][i] * cos(E->X[1][i]);
      }

 return;
 }

/*      ----------------------------------- */
  void  frames_for_dx(E,file_number)
  struct All_variables *E;
  int file_number; 
  {

  int i;
  static int nframe=0;
  FILE *fp;
  char output_file[255];
  const float offset1 = -1.2;
  const float offset2 = 0.2;

   nframe ++;

   sprintf(output_file,"%s/mv.%03d.dx",E->control.data_file,nframe);
   fp=fopen(output_file,"w");
   for (i=1;i<=E->mesh.nno;i++)
     fprintf(fp,"%g %g %g\n",E->ibm_dx.x1[i]+offset1,E->ibm_dx.x2[i],E->T[i]);
   fclose(fp);

   sprintf(output_file,"%s/nv.%03d.dx",E->control.data_file,nframe);
   fp=fopen(output_file,"w");
   for (i=1;i<=E->mesh.nno;i++)
     fprintf(fp,"%g %g %g\n",E->ibm_dx.x1[i]+offset2,E->ibm_dx.x2[i],E->C[i]);
   fclose(fp);

   sprintf(output_file,"%s/mv.%03d.general",E->control.data_file,nframe);
   fp=fopen(output_file,"w");
   fprintf(fp,"file = mv.%03d.dx\n",nframe);
   fprintf(fp,"grid = %2d x %2d\n",E->ibm_dx.nox,E->ibm_dx.noz);
   fprintf(fp,"format = ascii\n");
   fprintf(fp,"interleaving = field\n");
   fprintf(fp,"majority = row\n");
   fprintf(fp,"field = locations, field0\n");
   fprintf(fp,"structure = 2-vector, scalar\n");
   fprintf(fp,"type = float, float\n");
   fprintf(fp,"\n");
   fprintf(fp,"end\n");
   fclose(fp);

   sprintf(output_file,"%s/nv.%03d.general",E->control.data_file,nframe);
   fp=fopen(output_file,"w");
   fprintf(fp,"file = nv.%03d.dx\n",nframe);
   fprintf(fp,"grid = %2d x %2d\n",E->ibm_dx.nox,E->ibm_dx.noz);
   fprintf(fp,"format = ascii\n");
   fprintf(fp,"interleaving = field\n");
   fprintf(fp,"majority = row\n");
   fprintf(fp,"field = locations, field0\n");
   fprintf(fp,"structure = 2-vector, scalar\n");
   fprintf(fp,"type = float, float\n");
   fprintf(fp,"\n");
   fprintf(fp,"end\n");
   fclose(fp);



   return;
   }

void output_velo_related_binary(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int el,els,i,j,k,ii,m,node,fd;
  int nox,noz,noy,nfx,nfz,nfy1,nfy2,size1,size2;
  char output_file[255];
  static float *SV,*EV;
  static int been_here=0;

  void get_surface_velo ();
  void get_ele_visc ();
  const int nno = E->mesh.nno;
/*
  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.velo",E->control.data_file);
    E->filed[10]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_t",E->control.data_file);
    E->filed[11]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_b",E->control.data_file);
    E->filed[12]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.visc",E->control.data_file);
    E->filed[13]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.fas670",E->control.data_file);
    E->filed[14]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.stress",E->control.data_file);
    E->filed[9]=open(output_file,O_RDWR | O_CREAT, 0644);
    }

  if (been_here==0)  {
    ii = E->mesh.nsf;
    SV = (float *) malloc ((2*ii+2)*sizeof(float));

    size2 = (E->mesh.nel+1)*sizeof(float);
    EV = (float *) malloc (size2);
    been_here++;
    }

  ii = E->mesh.nsf;
  size2 = 2*(ii+2)*sizeof(float);
    get_surface_velo (E,SV);
  write(E->filed[10],SV,size2);

  size2 = (E->mesh.nsf+1)*sizeof(float);
  write(E->filed[11],E->slice.tpg,size2);
  write(E->filed[12],E->slice.tpgb,size2);

  size2 = (E->mesh.nel+1)*sizeof(float);
    get_ele_visc (E,EV);
  write(E->filed[13],EV,size2);

  size2 = (E->mesh.nsf+1)*sizeof(float);
  write(E->filed[14],E->Fas670_b,size2);

  size2 = (2*E->mesh.nsf+1)*sizeof(float);
  write(E->filed[9],E->stress,size2);
*/
  return;
  }

/* ====================================================================== */

void output_temp(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int nno,i,j,fd;
  static int *temp1;
  static int been_here=0;
  static int size2,size1;
  char output_file[255];
/*
  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.temp",E->control.data_file);
    E->filed[5]=open(output_file,O_RDWR | O_CREAT, 0644);
    }

  if (been_here==0) {
    temp1 = (int *) malloc ((E->mesh.noy*6)*sizeof(int));

    sprintf(output_file,"%s.mesh",E->control.data_file);
    E->filed[1]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.x",E->control.data_file);
    E->filed[2]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.z",E->control.data_file);
    E->filed[3]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.y",E->control.data_file);
    E->filed[4]=open(output_file,O_RDWR | O_CREAT, 0644);

    size1 = (E->mesh.noy*6)*sizeof(int);
    size2= (E->mesh.nno+1)*sizeof(float);

    temp1[1] = E->mesh.nno;
    temp1[3] = size2;
    temp1[5] = E->mesh.nsf;
    temp1[6] = E->mesh.nel;

        write(E->filed[1],temp1,size1);
        write(E->filed[2],E->X[1],size2);
        write(E->filed[3],E->X[2],size2);
        write(E->filed[4],E->X[3],size2);

    close(E->filed[1]);
    close(E->filed[2]);
    close(E->filed[3]);
    close(E->filed[4]);

    been_here++;
    }


    write(E->filed[5],E->T,size2);
*/
  return; 
}
