#include <stdio.h>
#include <math.h>
//#include <malloc.h>
#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void phase_change(E,Bb,Bb_b,Bt,Bt_b)
  struct All_variables *E;
  float *Bb,*Bb_b,*Bt,*Bt_b;
{
  static int been_here=0;

  FILE *fp1,*fp2;
  char output_file[255];
  void return_horiz_ave();

  int i,j,k,n,ns;
  double e_pressure,temp,temp1;

  static float pt5=0.5;
  static float one=1.0;

  if (been_here++==0)    {

      E->control.Ra_670 = E->control.Ra_670*E->control.Ra_temp
              /(E->data.density*E->data.therm_exp*E->data.ref_temperature);
      E->control.Ra_410 = E->control.Ra_410*E->control.Ra_temp
              /(E->data.density*E->data.therm_exp*E->data.ref_temperature);

      E->control.clapeyron670 = E->control.clapeyron670*E->data.ref_temperature/
                          (E->data.density*E->data.grav_acc*E->sphere.ro_dim);
      E->control.clapeyron410 = E->control.clapeyron410*E->data.ref_temperature/
                          (E->data.density*E->data.grav_acc*E->sphere.ro_dim);

      E->control.width670 = E->sphere.ro_dim/E->control.width670;
      E->control.width410 = E->sphere.ro_dim/E->control.width410;

      E->control.transT670 = E->control.transT670/E->data.ref_temperature;
      E->control.transT410 = E->control.transT410/E->data.ref_temperature;

fprintf(E->fp,"Rab410 670=%g %g Clap410 670=%g %g %g %g %g %g\n",E->control.Ra_410,E->control.Ra_670,E->control.clapeyron410,E->control.clapeyron670,E->control.width670,E->control.width410,E->control.transT670, E->control.transT410);
fflush (E->fp);
      }

  return_horiz_ave(E,E->T,E->Have.T);

  temp1 = 0.0;
  for(i=1;i<E->mesh.noz;i++)  {
     if (E->viscosity.zlm<=E->X[2][i+1]&&E->viscosity.zlm>=E->X[2][i])  {
         temp1 = E->Have.T[i] + (E->Have.T[i+1]-E->Have.T[i])*(E->viscosity.zlm-E->X[2][i])/(E->X[2][i+1]-E->X[2][i]);
         break;
         }
     }
  E->control.transT670 = temp1;
  temp1 = 0.0;
  for(i=1;i<E->mesh.noz;i++)  {
     if (E->viscosity.z410<=E->X[2][i+1]&&E->viscosity.z410>=E->X[2][i])  {
         temp1 = E->Have.T[i] + (E->Have.T[i+1]-E->Have.T[i])*(E->viscosity.z410-E->X[2][i])/(E->X[2][i+1]-E->X[2][i]);
         break;
         }
     }
  E->control.transT410 = temp1;


  for(i=1;i<=E->mesh.nno;i++)  {
    e_pressure = E->viscosity.zlm - E->X[2][i] -
            E->control.clapeyron670*(E->T[i]-E->control.transT670);
    Bb[i] = pt5*(one+tanh(E->control.width670*e_pressure));
    }

  for(i=1;i<=E->mesh.nno;i++)  {
    e_pressure = E->viscosity.z410 - E->X[2][i] -
            E->control.clapeyron410*(E->T[i]-E->control.transT410);
    Bt[i] = pt5*(one+tanh(E->control.width410*e_pressure));
    }


if (E->advection.timesteps%E->control.record_every == 0)   {

    for (j=1;j<=E->mesh.nox;j++)  {
      Bb_b[j]=0.0;
      for (i=1;i<=E->mesh.noz;i++)   {
        n = (j-1)*E->mesh.noz + i;
        if (Bb[n]>=pt5&&Bb[n+1]<=pt5)   {
          Bb_b[j]=(E->X[2][n+1]-E->X[2][n])*(pt5-Bb[n])/(Bb[n+1]-Bb[n])+E->X[2][n];
          break;
          }
        }
      }

   for (j=1;j<=E->mesh.nox;j++)  {
      Bt_b[j]=0.0;
      for (i=1;i<=E->mesh.noz;i++)   {
        n = (j-1)*E->mesh.noz + i;
        if (Bt[n]>=pt5&&Bt[n+1]<=pt5)  {
          Bt_b[j]=(E->X[2][n+1]-E->X[2][n])*(pt5-Bt[n])/(Bt[n+1]-Bt[n])+E->X[2][n];
          break;
          }
        }
      }



  sprintf(output_file,"%s/fas.%d",E->control.data_file,E->advection.timesteps);
  fp1=fopen(output_file,"w");
  for (j=1;j<=E->mesh.nox;j++)
    fprintf(fp1,"%.4e %.5e %.5e\n",E->X[1][j*E->mesh.noz],Bt_b[j],Bb_b[j]);
  fclose(fp1);

   }

if (E->monitor.solution_cycles%E->control.record_every == 0)   {
    fprintf(E->fp,"fas=%g %g %g %g %g %g %g %g\n", E->control.clapeyron410,E->control.clapeyron670 ,E->control.Ra_410,E->control.Ra_670,E->control.transT410,E->control.transT670,E->control.width410,E->control.width670);fflush(E->fp);
	        }



  return;
  }


/********************************************************
 * melting_by_node                                      *
 *			funciton to compute and track the mantle				*
 *			melting when the temperature exceeds the				*
 *			solidus.	Melt fraction is determined for each	*
 *			element.  Melt is assumed to instantly be				*
 *			removed to the surface.  Composition of element	*
 *			based on ratio of pristine mantle to residuum.	*
 *			Compositional tracers re-allocated based on			*
 *			new composition.																*
 *                                                      *
 *			Will need to deal with the situation in which		*
 *			element is entirely melted.  Also compaction of	*
 *			pore space when melt removed.  Propogate solids	*
 *			down.  Topo is solid/liquid interface?					*
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *			dT			Temp increase, assuming no latent heat	*
 *			n				node																		*
 *			k				radial node															*
 ********************************************************/

float melting_by_node(E,dT,n,k)
  struct All_variables *E;
	float dT;
	int n;
	int k;
{
	char output_file[50];

	float Ed;			/* Energy density */
	float Hf;			/* Latent heat of fusion */
	float beta;		/* dF/dT */
	float dTsol;	/* dT in excess of solidus */
  float Tsl;    /* difference btw liquidus and solidus */
	float F;			/* Melt fraction */
	float dF;			/* Change in melt fraction, assuming C = 0 */
	float Tboil;	/* Boiling point */
	int meltcase;	

	Hf = 4.187e5;			/* J/kg */
	beta = 0.00093;		/* K-1 */
	Tboil = 3223.0;		/* K, SiO2 boils, MgO still liquid */
  Tsl = E->liquidus[k] - E->solidus[k];

	/*	Convert dT to equivalent energy density */
	Ed = E->data.Cp*dT;
	dTsol = E->T[n] + dT - E->solidus[k];
//	fprintf(E->fp,"%d %f ",n,dTsol);

	/* Nondim stuff */
	Hf /= E->data.DeltaT;
	Tboil = Tboil/E->data.DeltaT - E->data.surf_temp;
	beta *= E->data.DeltaT;

	/* Three possible cases:
	 * 1. Tn below solidus, no melting, Tn = To + dT
	 * 2. Tn above liquidus, complete melting, Tn = To + dT - L
	 * 3. Tn between liquidus and solidus, partial melting.
	 * Parameterization from Katz et al. (2003), G3 4, 1073.
	 * For now, I'm assuming a linear relationship between Temperature
	 * and melt fraction.
	 */

	if( dTsol < 0.0 ) {	/* Below Solidus */
			E->T[n] += dT;
			F = 0.0;
			meltcase=0;
		} 
		//else if ( (E->T[n]+dT-(Hf/E->data.Cp)) > E->lherzliq[k] ) {
		//else if ( dTsol > (1.0/beta + Hf/E->data.Cp) ) { /* WRONG */
		else if ( dTsol > (Tsl + Hf/E->data.Cp) ) { /* Above Liquidus */
			E->T[n] += dT - Hf/E->data.Cp;
			F = 1.0 - E->C[n];
			meltcase=1;

			/* T cannot rise above boiling point.  Ignore vapor. */
//			if(E->T[n] > Tboil)
//				E->T[n] = Tboil;
			if(E->T[n] > E->liquidus[k])
				E->T[n] = E->liquidus[k];
		} 
		else {	/* Partial Melt */
			/*
			E->T[n] += ( Ed*(E->lherzliq[k]-E->solidus[k]) 
																			- latent*(E->T[n]-E->solidus[k]) ) 
									/ ( latent + E->data.Cp*(E->lherzliq[k]-E->solidus[k]) );
			F = (E->T[n] - E->solidus[k]) / (E->lherzliq[k] - E->solidus[k]);*/
		
    /* Method of Elkins-Tanton (2005) -- DO NOT USE */
		/*	dF = dTsol / (Hf/E->data.Cp + 1.0/beta);*/		/* ratio of energy to latent heat */
     
     dF = (dTsol + E->C[n]*Hf/E->data.Cp) / (Tsl + Hf/E->data.Cp);  /* Frac of melt */
     if (dF > 1.0) dF = 1.0;  /* Can't have more than 100% melting */

		 if (dF > E->C[n])
				F = (dF - E->C[n]);   /* New melt created */
				/*F = (dF - E->C[n]) / (1 - E->C[n]); */
		 else 
				F = 0.0;
        
			meltcase=2;

/*			E->T[n] += dT - dTsol + F/beta;*/
      E->T[n] += dT - F*Hf/E->data.Cp;

		}

  if(E->C[n]+F > 1.0)
		fprintf(E->fp,"meltcase %d %d %f %f\n",n,meltcase,F,E->C[n]);
//		E->Fm[n] = F;  /* This is returned, don't do it here!

	/* Determine composition of element from melt fraction
		 Element can never become less melted */
/*		if (F > E->C[n])  *//* unnecessary? */
			E->C[n] += F;  /* effectively, C = dF */

//		fprintf(E->fp,"\n");
		return(F);
}


/********************************************************
 * melting				                                      *
 *			funciton to compute and track the mantle				*
 *			melting when the temperature exceeds the				*
 *			solidus.	Melt fraction is determined for each	*
 *			element.  Melt is assumed to instantly be				*
 *			removed to the surface.  Composition of element	*
 *			based on ratio of pristine mantle to residuum.	*
 *			Compositional tracers re-allocated based on			*
 *			new composition.																*
 *                                                      *
 *			Will need to deal with the situation in which		*
 *			element is entirely melted.  Also compaction of	*
 *			pore space when melt removed.  Propogate solids	*
 *			down.  Topo is solid/liquid interface?					*
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *			dT			Temp increase, assuming no latent heat	*
 *			e     	element																	*
 *			k				radial element													*
 ********************************************************/

float melting(E,dT,e,k)
  struct All_variables *E;
	float dT;
	int e;
	int k;
{
	char output_file[50];

	float Ed;			/* Energy density */
	float Hf;			/* Latent heat of fusion */
	float beta;		/* dF/dT */
	float dTsol;	/* dT in excess of solidus */
  float Tsl;    /* difference btw liquidus and solidus */
	float F;			/* Melt fraction */
	float dF;			/* Change in melt fraction, assuming C = 0 */
	float Tboil;	/* Boiling point */
  float elemsol, elemliq;  /* Solidus and Liquidus of element */
	int meltcase;	

	Hf = 4.187e5;			/* J/kg */
	beta = 0.00093;		/* K-1 */
	Tboil = 3223.0;		/* K, SiO2 boils, MgO still liquid */

  elemsol = (E->solidus[k] + E->solidus[k+1]) / 2.0;
  elemliq = (E->liquidus[k] + E->liquidus[k+1]) / 2.0;
  Tsl = elemliq - elemsol;

	/*	Convert dT to equivalent energy density */
	Ed = E->data.Cp*dT;
	dTsol = E->TE[e] + dT - elemsol;
//	fprintf(E->fp,"%d %f ",n,dTsol);

	/* Nondim stuff */
	Hf /= E->data.DeltaT;
	Tboil = Tboil/E->data.DeltaT - E->data.surf_temp;
	beta *= E->data.DeltaT;

	/* Three possible cases:
	 * 1. Tn below solidus, no melting, Tn = To + dT
	 * 2. Tn above liquidus, complete melting, Tn = To + dT - L
	 * 3. Tn between liquidus and solidus, partial melting.
	 * Parameterization from Katz et al. (2003), G3 4, 1073.
	 * For now, I'm assuming a linear relationship between Temperature
	 * and melt fraction.
	 */

	if( dTsol < 0.0 ) {	/* Below Solidus */
			E->TE[e] += dT;
			F = 0.0;
			meltcase=0;
		} 
		else if ( dTsol > (Tsl + Hf/E->data.Cp) ) { /* Above Liquidus */
		/*else if ( dTsol > (1.0/beta + Hf/E->data.Cp) ) {*/ /* WRONG */
			E->TE[e] += dT - Hf/E->data.Cp;
			F = 1.0 - E->CE[e];
			meltcase=1;

			/* T cannot rise above boiling point.  Ignore vapor. */
			if(E->TE[e] > elemliq)
				E->TE[e] = elemliq;
		} 
		else {	/* Partial Melt */
			/* Method of Elkins-Tanton (2005) DO NOT USE */
/*			dF = dTsol / (Hf/E->data.Cp + 1.0/beta);*/		/* ratio of energy to latent heat */
/*			if (dF > E->CE[e])
				F = (dF - E->CE[e]); */ /* New melt created */
/*			else 
				F = 0.0;*/
			meltcase=2;

      dF = (dTsol + E->CE[e]*Hf/E->data.Cp) / (Tsl + Hf/E->data.Cp);  /* Frac of melt */
     
     if (dF > 1.0) dF = 1.0;  /* Can't have more than 100% melting */
		 if (dF > E->CE[e])
				F = (dF - E->CE[e]);   /* New melt created */
		 else 
				F = 0.0;
        
			/*E->TE[e] += dT - dTsol + F/beta;*/
      E->TE[e] += dT - F*Hf/E->data.Cp;
		}

//		fprintf(E->fp,"%d %f ",meltcase,F);

	/* Determine composition of element from melt fraction
		 Element can never become less melted */
/*		if (F > E->C[n])  *//* unnecessary? */
			E->CE[e] += F;  /* effectively, C = dF */

//		fprintf(E->fp,"\n");
		return(F);
}


/********************************************************
 * freezemelt                                           *
 *      function to compute the amount of freezing of   *
 *      the subsurface ocean or melting at the base of  *
 *      the ice shell based on the mismatch between the *
 *      heat produced by the core over the current      *
 *      timestep and the heat flux across the base of   *
 *      the ice shell.                                  *
 *                                                      *
 * Parameters                                           *
 *      E       All_variables                           *
 *      Fb      Basal heat flux                         *
 *                                                      *
 * Returns                                              *
 *      deltarb Change in ice shell thickness           *
 ********************************************************/

float freezemelt(E, Fb)
  struct All_variables *E;
  float Fb;
{

  float dimfactor;  /* Prefactor relating thickness change to heat flux */
  float Fc;         /* Heat flux coming out of the core at base of ice shell */
  float deltarb;    /* Change in shell thickness */

  dimfactor = (E->data.Cp * E->data.DeltaT) / E->data.Hf;

  Fc = (E->control.Qc/3.0) * pow(E->sphere.rcore,3.0) / pow(E->sphere.ri,2.0);

  deltarb = dimfactor * (Fb - Fc) * E->advection.timestep;

  /* Limit deltarb based on space available. Can't freeze past core */
  deltarb = min(deltarb,(E->sphere.ri-E->sphere.rcore));

  fprintf(E->fp,"%e Fb %e Fc %e dt %e\n",dimfactor,Fb,Fc,E->advection.timestep);
  fprintf(E->fp,"Change in shell thickness: %e = %e m\n",deltarb,
                                    (deltarb*E->sphere.ro_dim));

//  fprintf(stderr,"%e\n",deltarb);
  return deltarb;
}
