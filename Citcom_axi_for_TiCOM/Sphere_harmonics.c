#include <math.h>
#include <fcntl.h>
#include <sys/types.h>
#include "global_defs.h"
#include "element_definitions.h"

const double shl[5][5]; 
                  

/* ====================================================*/
 void sphere_harmonics_initialinaztion(E)
     struct All_variables *E;
 {
    int i,ll,mm,lll;
 void set_sphereh_consts_index();
 void  compute_sphereh_table();


  input_int ("llmax",&(E->sphere.llmax),"10");
  input_int ("selx",&(E->sphere.elx),"10");

    E->sphere.nox = E->sphere.elx + 1;
    E->sphere.nsf = E->sphere.nox;
    E->sphere.snel = E->sphere.elx;

  E->sphere.output_llmax = E->sphere.llmax;

  E->sphere.sx[1] = (double *) malloc((E->sphere.nsf+1)*sizeof(double));
  E->sphere.sx[2] = (double *) malloc((E->sphere.nsf+1)*sizeof(double));
  for (i=0;i<=E->sphere.llmax;i++)
    E->sphere.hindex[i] = i;

  E->sphere.hindice = E->sphere.llmax;

  E->sphere.con  = (double *)malloc((E->sphere.hindice+3)*sizeof(double));
  E->sphere.sphc = (double *)malloc((E->sphere.hindice+3)*sizeof(double));

  E->sphere.field = (double *)malloc((E->sphere.nsf+3)*sizeof(double));

  for (i=1;i<=E->sphere.elx*2;i++)
    E->sphere.tableplm[i]= (double *)malloc((E->sphere.hindice+3)*sizeof(double));

  for (i=1;i<=E->sphere.nox;i++)
    E->sphere.tableplm_n[i]= (double *)malloc((E->sphere.hindice+3)*sizeof(double));

  E->sphere.sien  = (struct SIEN *) malloc((E->sphere.snel+1)*sizeof(struct SIEN));

  set_sphereh_consts_index(E);

  compute_sphereh_table(E);

 return;
 }


/*  ======================================================================
    ======================================================================  */
 void set_sphereh_consts_index(E)
   struct All_variables *E;

  {
  int rr,ends,es,node,ll,lll,mm,i,j;
  double dth,dfi,sqrt_multis(),gauss;


  ends = 2;

  for (i=1;i<=E->sphere.elx;i++) {
    es = i;
    node = i;
    for (rr=1;rr<=ends;rr++)
        E->sphere.sien[es].node[rr] = node 
		       + offset[rr].vector[1];
    }

  for (ll=0;ll<=E->sphere.output_llmax;ll++)   {
     E->sphere.con[E->sphere.hindex[ll]] = 
	     sqrt( (2*ll+1)/(4.0*M_PI) )
	    *sqrt_multis(ll,ll);  /* which is sqrt((ll-mm)!/(ll+mm)!) */
     }

  dth = M_PI/E->sphere.elx;

  for (i=1;i<=E->sphere.nox;i++) {
    E->sphere.sx[1][i] = dth*(i-1);
    }

  gauss = 1.0/sqrt(3.0);

/* for a=1, intp=1-2 */
  shl[1][1] = 0.5*(1+gauss);
  shl[1][2] = 0.5*(1-gauss);
/* for a=2, intp=1-2 */
  shl[2][1] = 0.5*(1-gauss);
  shl[2][2] = 0.5*(1+gauss);

  return;
  }
  

/* ==================================================*/
 void  compute_sphereh_table(E)
 struct All_variables *E;
 {

 int lll,intp,rr,node,ends,ll,mm,es,i,j,p,jj;
 double t,f,modified_plgndr_a();
 const double pt25=0.25;


  for (i=1;i<=E->sphere.elx;i++) {
    es = i;
    for (intp=1;intp<=2;intp++)   {
      t = shl[1][intp]*E->sphere.sx[1][E->sphere.sien[es].node[1]]
        + shl[2][intp]*E->sphere.sx[1][E->sphere.sien[es].node[2]];
      jj = (i-1)*2+intp;
      E->sphere.tablesint[jj] = sin(t);
      for (ll=0;ll<=E->sphere.output_llmax;ll++)   {
         mm = 0;
         p = E->sphere.hindex[ll];
         E->sphere.tableplm[jj][p] = modified_plgndr_a(ll,mm,t) ;
         }
      } 
    } 

  for (i=1;i<=E->sphere.nox;i++) {
    node = i;
    t=E->sphere.sx[1][node];
    E->sphere.tablesint_n[i] = sin(t);
    for (ll=0;ll<=E->sphere.output_llmax;ll++)  {
         mm=0;
         p = E->sphere.hindex[ll];
         E->sphere.tableplm_n[i][p] = modified_plgndr_a(ll,mm,t) ;
         }
    }


 return;
 }

/* =========================================================
  ========================================================= */
 void raw_sphere_expansion(E,X1,X2,ii,iprint)
 struct All_variables *E;
 float *X1, *X2;
 int iprint,ii;
 {
 int intp,p,i,j,es,mm,ll,lll,i1,i2,i3,i4;
 double temp,temp1,area,t,f,*TG,*sphc;
 const double pt25=0.25;
 static int been_here=0;
 void sphere_expansion();
 void output_field_spectral();
 void inv_sphere_harmonics();
 FILE *fp;


   TG = (double *)malloc((E->sphere.nox+1)*sizeof(double));
   sphc= (double *)malloc((E->sphere.llmax+3)*sizeof(double));

   for (i=0;i<=E->sphere.hindice;i++)    {
      sphc[i] = 0.0;
      }

   i1  = 1;
   TG[E->sphere.nox] = X2[E->mesh.nox];

   for (i=1;i<E->sphere.nox;i++)            {
     for (j=i1;j<=E->mesh.elx;j++)  {
       if (E->sphere.sx[1][i]>=X1[j] && E->sphere.sx[1][i]<=X1[j+1])    {
          TG[i] = (X2[j+1]-X2[j])/(X1[j+1]-X1[j])*(E->sphere.sx[1][i]-X1[j]) + X2[j];
          i1 = j;
          break;
          }
       }
   }

  sphere_expansion(E,TG,sphc);

  if (iprint)
     output_field_spectral(E,sphc,TG,ii+1,0);
    
  free ((void *)TG);
  free ((void *)sphc);

  return;
  }

/* =========================================================
  ========================================================= */
 void sphere_expansion(E,TG,sphc)
 struct All_variables *E;
 double *TG,*sphc;
 {
 int intp,p,i,j,es,mm,ll,lll,i1,i2,i3,i4;
 double temp,temp1,area,t,f,*TGAG[5];
 const double pt25=0.25;
 static int been_here=0;
 FILE *fp;

   for (intp=1;intp<=2;intp++)
     TGAG[intp] = (double *)malloc((E->sphere.elx+1)*sizeof(double));


   for (i=0;i<=E->sphere.hindice;i++)    {
      sphc[i] = 0.0;
      }

   for (i=1;i<=E->sphere.elx;i++) {

        es = i;
        i1 = E->sphere.sien[es].node[1];
        i2 = E->sphere.sien[es].node[2];

        for (intp=1;intp<=2;intp++)
           TGAG[intp][es] = TG[i1]*shl[1][intp]+TG[i2]*shl[2][intp];
        }


   area = 0.5*2.0*M_PI*M_PI/(E->sphere.elx);

   for (ll=0;ll<=E->sphere.output_llmax;ll++)   {

     mm = 0;
     p = E->sphere.hindex[ll];

     for (i=1;i<=E->sphere.elx;i++) {
	es = i;
        for (intp=1;intp<=2;intp++)  
             sphc[p]+=TGAG[intp][es]*E->sphere.tableplm[(i-1)*2+intp][p]*E->sphere.tablesint[(i-1)*2+intp]; 
	 }  

     sphc[p] *= area; 

     }       /* end for ll and mm  */


   for (intp=1;intp<=2;intp++)
     free ((void *)TGAG[intp]);

 return;
 }

 /* =========================================================== */
 void inv_sphere_harmonics(E,sphc,TG)
 struct All_variables *E;
 double *TG,*sphc;
 {
 int lll,k,ll,mm,node,i,j,p,noz,snode;

   for (i=1;i<=E->sphere.nox;i++)  
       TG[i]=0.0;

   for (ll=0;ll<=E->sphere.output_llmax;ll++)   {
     p = E->sphere.hindex[ll];

     for (node=1;node<=E->sphere.nox;node++)    
       TG[node] += (sphc[p]*E->sphere.tableplm_n[node][p]);
     }

 return;
 }

/* ====================================================*/
/* ====================================================*/
void output_field_spectral(E,sphc,TG,ii,icon)
     struct All_variables *E;
     double *sphc,*TG;
     int icon,ii;
{
   FILE * fp;
   char filename[100];
   int lll,i,j,ll,mm,node,noy;
   double tmin,tmax,t,f,rad;
   double *TG1;

   rad = 180.0/M_PI;

   sprintf(filename,"%s.spec.%d",E->control.data_file,ii);
   fp = fopen(filename,"w");

   TG1 = (double *)malloc((E->sphere.nox+1)*sizeof(double));

     t = 0.0;
     for (node=1;node<E->sphere.nox;node++)
         t += TG[node] + TG[node+1];

     t = 0.5*t/E->sphere.elx;

     tmax = -1e10;
     tmin = 1e10;
     for (node=1;node<=E->sphere.nox;node++)    {
         TG1[node] = TG[node] - t;
         if (TG1[node]>tmax)tmax = TG1[node];
         if (TG1[node]<tmin)tmin = TG1[node];
         }

     fprintf(fp,"%.7e\n",tmax-tmin);
     for (node=1;node<=E->sphere.nox;node++)
         fprintf(fp,"%3d %.7e %.7e\n",node,TG[node],TG1[node]);

     for (ll=0;ll<=E->sphere.output_llmax;ll++)   {
         i = E->sphere.hindex[ll];
         fprintf(fp,"%3d %.12e \n",ll,sphc[i]);
         }

   fclose(fp);

  return;
  }

/* =========================================================
  ========================================================= */
 void test(E,X1,X2,ii,iprint)
 struct All_variables *E;
 float *X1, *X2;
 int iprint,ii;
 {
 int intp,p,i,j,es,mm,ll,lll,i1,i2,i3,i4;
 double temp,temp1,area,t,f,*TG,*sphc;
 const double pt25=0.25;
 static int been_here=0;
 void sphere_expansion();
 void output_field_spectral();
 void inv_sphere_harmonics();
 FILE *fp;

 fp = fopen("junk","w");

   TG = (double *)malloc((E->sphere.nox+1)*sizeof(double));
   sphc= (double *)malloc((E->sphere.llmax+3)*sizeof(double));

   for (i=0;i<=E->sphere.hindice;i++)    {
      sphc[i] = 0.0;
      }

  ll= E->convection.perturb_ll[0];
  sphc[ll] = E->convection.perturb_mag[0];

  inv_sphere_harmonics(E,sphc,TG);

   for (i=1;i<E->sphere.nox;i++)
      fprintf(fp,"%d %.12e\n",i,TG[i]);

  sphere_expansion(E,TG,sphc);

  if (iprint)
     output_field_spectral(E,sphc,TG,ii,0);
    


   i1  = 1;
   TG[E->sphere.nox] = X2[E->mesh.nox];

   for (i=1;i<E->sphere.nox;i++)            {
     for (j=i1;j<=E->mesh.elx;j++)  {
       if (E->sphere.sx[1][i]>=X1[j] && E->sphere.sx[1][i]<=X1[j+1])    {
          TG[i] = (X2[j+1]-X2[j])/(X1[j+1]-X1[j])*(E->sphere.sx[1][i]-X1[j]) + X2[j];
          i1 = j;
          break;
          }
       }
   }

   for (i=1;i<E->sphere.nox;i++)
      fprintf(fp,"%d %.12e\n",i,TG[i]);

  sphere_expansion(E,TG,sphc);

  if (iprint)
     output_field_spectral(E,sphc,ii+1,0);
    
  free ((void *)TG);
  free ((void *)sphc);

  return;
  }
