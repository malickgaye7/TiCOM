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
#include <complex.h>

#ifdef _UNICOS
#include <fortran.h>
#endif

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

int epsilon[4][4] = {   /* Levi-Cita epsilon */
  {0, 0, 0, 0},
  {0, 1,-1, 1},
  {0,-1, 1,-1},
  {0, 1,-1, 1} };

static float cost_per_level[MAX_LEVELS]; /* this will accumulate data over the run */
static int total_cycles[MAX_LEVELS];


/*=====================================================================
  Variable dimension matrix allocation  function from numerical recipes
  Note: ANSII consistency requires some additional features !
  =====================================================================  */
 
double **dmatrix(nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch;
{
  int i,nrow = nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointer to rows  */
  m=(double **) malloc((nrow+1)* sizeof(double *));
  m+=1;
  m-= nrl;

  /*  allocate rows and set the pointers accordingly   */
  m[nrl] = (double *) malloc((nrow*ncol+1)* sizeof(double));
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)
     m[i] = m[i-1] + ncol; 

  return(m);		}


float **fmatrix(nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch;
{
  int i,nrow = nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointer to rows  */
  m=(float **) malloc((unsigned)((nrow+1)* sizeof(float *)));
  m+=1;
  m-= nrl;

  /*  allocate rows and set the pointers accordingly   */
  m[nrl] = (float *) malloc((unsigned)((nrow*ncol+1)* sizeof(float)));
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)
     m[i] = m[i-1] + ncol; 

  return(m);		}


void dfree_matrix(m,nrl,nrh,ncl,nch)
     double **m;
     int nrl,nrh,ncl,nch;
{	
  int i;
  for(i=nrh;i>=nrl;i--) 
    free((void *)(m[i] + ncl));
  free((void *) (m+nrl));
  return;	
}

void ffree_matrix(m,nrl,nrh,ncl,nch)
     float **m;
     int nrl,nrh,ncl,nch;
{	
  int i;
  for(i=nrh;i>=nrl;i--) 
    free((void *)(m[i] + ncl));
  free((void *) (m+nrl));
  return;	
}

/*=============================================================
  Functions to allocate/remove space for variable sized vector.
  =============================================================  */

double *dvector(nl,nh)
     int nl,nh;
{	
  double *v;
  v=(double *) malloc((unsigned) ( nh - nl +1)* sizeof(double));
  return( v-nl );  }

float *fvector(nl,nh)
     int nl,nh;
{	
  float *v;
  v=(float *) malloc((unsigned) ( nh - nl +1)* sizeof(float));
  return( v-nl );  }

void dfree_vector(v,nl,nh)
     double *v;
     int nl,nh;
{	
  free((char*) (v+nl));	}

void ffree_vector(v,nl,nh)
     float *v;
     int nl,nh;
{	
  free((char*) (v+nl));	}

int *sivector(nl,nh)
     int nl,nh;
{
  int *v;
  v=(int*) malloc((unsigned)(nh-nl +1) * sizeof(int));
  return (v-nl);
}

void sifree_vector(v,nl,nh)
     int *v;
     int nl,nh;
{ free((char *) (v+nl));    }


double pdot(E,A,B,lev)
     struct All_variables *E;
     double *A,*B;
     int lev;
    
{   double prod;
    int e;
    const int n=E->mesh.NPNO[lev];

    prod = 0.0;
    for(e=1;e<=n;e++)
      prod += A[e] * B[e] ;
    return(prod);  
    }

double pselfdot(E,A)
     struct All_variables *E;
     double *A;
    
{   double prod;
    int e;
    const int n=E->mesh.npno;

    prod = 0.0;
    for(e=1;e<=n;e++)
      prod += A[e] * A[e] ;
    return(prod);  }

/* Beware of the alias if A=B on vector machines, use vselfdot instead  */
float fvdot(E,A,B,level)
     struct All_variables *E;
     float *A,*B;
     int level;
{   
    float prod,mprod[1];
    int i,incx=1;
    
    char trans='N';
    double alpha=1.0;
    
    const int n = E->mesh.NEQ[level];

#ifdef blas
   
    
#else    
    prod = 0.0;
    for(i=0;i<n;i++)
	prod += A[i] * B[i] ;
#endif

    return(prod);  }

double vdot(E,A,B,level)
     struct All_variables *E;
     double *A,*B;
     int level;
{   
    double prod,mprod[1];
    int i,incx=1;
    
    char trans='N';
    double alpha=1.0;
    
    const int n = E->mesh.NEQ[level];

#ifdef blas
   
    
#else    
    prod = 0.0;
    for(i=0;i<n;i++)
	prod += A[i] * B[i] ;
#endif

    return(prod);  }


double vselfdot(E,A,level)  
     struct All_variables *E;
     double *A;
     int level;
{ 
    double prod;
    int i,n;

  n = E->mesh.NEQ[level];

  prod = 0.0;
  for(i=0;i<n;i++)
    prod += A[i] * A[i] ;
  
  return(prod);
}


double vfselfdot(E,A,level)  
     struct All_variables *E;
     float *A;
     int level;
    

{ float prod;
  int i,n;

     
  n = E->mesh.NEQ[level];

  prod = 0.0;
  for(i=0;i<n;i++)
    prod += A[i] * A[i] ;
  

  return(prod);  }

float fdot(A,B,n1,n2)
     float *A,*B;
     int n1,n2;

{  float prod;
   int i;
   
   prod = 0.0;
   for(i=n1;i<=n2;i++)
     prod += A[i] * B[i] ;
  

   return(prod);  }

float fselfdot(A,n1,n2)
     float *A;
     int n1,n2;

{  float prod;
   int i;
 
   prod = 0.0;
   for(i=n1;i<=n2;i++)
     prod += A[i] * A[i] ;

   return(prod);  }


float dot(E,A,B)
     struct All_variables *E;
     float *A,*B;

{   float prod = 0.0;
    float domega;
    int e,i,j;

    for(e=1;e<=E->mesh.nel;e++)
      for(i=1;i<=vpoints[E->mesh.nsd];i++)
	{ j=E->ien[e].node[i];
	  domega =  E->ECO[E->mesh.levmax][e].area;
	  prod += A[j] * B[j] * domega; }

    return(prod);  }

float selfdot(E,A)
     struct All_variables *E;
     float *A;

{   double prod = 0.0;
    double domega;
    int e,i,j;

    for(e=1;e<=E->mesh.nel;e++)
      for(i=1;i<=vpoints[E->mesh.nsd];i++)
	{ j=E->ien[e].node[i];
	  domega =  E->ECO[E->mesh.levmax][e].area;
	  prod += A[j] * A[j] * domega; }

    return(prod);  }

void dvcopy(A,B,a,b)
     double *A,*B;
     int a,b;

{   int i;

    for(i=a;i<=b;i++)
      A[i] = B[i]; 

    return; }

void vcopy(A,B,a,b)
     float *A,*B;
     int a,b;

{   int i;

    for(i=a;i<=b;i++)
      A[i] = B[i]; 

    return; }


void vprod(R,A,B,a,b)
     double *R,*A,*B;
     int a,b;

{   int i;

    for(i=a;i<=b;i++)
      R[i] = A[i] * B[i]; 

    return; }

float fnmax(E,A,a,b)
     struct All_variables *E;
     float *A;
     int a,b;

{ float maxm = -1.0e32;
  int i;

  for(i=a;i<=b;i++)
    if (A[i]>maxm) maxm = A[i];

  return(maxm);
}




/*  ===========================================================
    Iterative solver also using multigrid  ........
    ===========================================================  */

int solve_del2_u(E,d0,F,acc,high_lev)    
     struct All_variables *E;
     double *d0;
     double *F;
     double acc;
     int high_lev;
{ 
    void assemble_del2_u();
    void gauss_seidel();
    void e_assemble_del2_u();
    void n_assemble_del2_u();
    void strip_bcs_from_residual();
   
    double conj_grad();
    double multi_grid();
    double vdot(),vselfdot();
   
    static int been_here = 0;
    int count,cycles,convergent;
    int i, neq, gneq;
  
    double CPU_time0(),initial_time,time;
    double acc1,residual,prior_residual,r0;
    double *D1, *r, *Au;
  
    neq  = E->mesh.NEQ[high_lev];
    gneq  = E->mesh.NEQ[high_lev];
    r    = (double *) malloc((neq + 10)*sizeof(double));
    D1   = (double *) malloc((neq + 10)*sizeof(double)); 
    Au   = (double *) malloc((neq + 10)*sizeof(double));
  
    if (been_here==0)  {
	E->control.total_iteration_cycles = 0;
	E->control.total_v_solver_calls = 0;
	for(i=E->mesh.levmin;i<=E->mesh.levmax;i++) {
	    cost_per_level[i] = 0.0;
	    total_cycles[i] = 0; 
	    }
	been_here++;
        }
  
    for(i=0;i<neq;i++)  {
	r[i] = F[i];
	Au[i] = d0[i] = D1[i] = 0.0;
	}
    
    r0=residual=sqrt(vdot(E,r,r,high_lev)/gneq);

    acc = max(acc,r0*E->control.accuracy);
   
    prior_residual=2*residual; 
    count = 0;  
    initial_time=CPU_time0();

    if (!(E->control.NMULTIGRID || E->control.EMULTIGRID)) {
	    cycles = E->control.v_steps_low;
	    time=CPU_time0();
	    residual = conj_grad(E,D1,r,Au,acc,&cycles,high_lev); 
	    for(i=0;i<neq;i++)  {
		d0[i] += D1[i];
		D1[i] = 0.0;
	    }

/*        count =0;
        do {
            gauss_seidel(E,D1,r,Au,acc,&cycles,high_lev,0);
            for(i=0;i<neq;i++) {
                r[i] = r[i]-Au[i];
                d0[i] += D1[i];
                }
            count ++;
            residual=sqrt(vdot(E,r,r,high_lev)/gneq);
            fprintf(stderr,"resi = %.6e for iter %d \n",residual,count);
            } while (residual > acc);
*/
	    cost_per_level[high_lev] +=  CPU_time0()-time;
	    total_cycles[high_lev] += cycles;
	 
	}
      
   else  { /* multigrid  */

        count =0;
        fprintf(E->fp,"resi = %.6e for iter %d acc %.6e\n",residual,count,acc);
        //fflush(E->fp); // DEBUG
        do {

	    residual=multi_grid(E,D1,r,Au,acc,high_lev);
	  
	      for(i=0;i<neq;i++)   {
		d0[i] += D1[i];
		D1[i] = 0.0;
       	        }
            count ++;
/*            fprintf(E->fp,"resi = %.6e for iter %d acc %.6e\n",residual,count,acc); */
/*            } while (residual > acc && count<E->control.max_vel_iterations);
 */   
            } while (residual > acc);
	}
	

	if((count > 0) && 
	   (residual > r0*2.0)  ||
	   (fabs(residual-prior_residual) < acc*0.1 && (residual > acc * 10.0))   )
	    convergent=0; 
	else {
	    convergent=1;
	    prior_residual=residual;
	    }
	
    if(E->control.print_convergence)   {
	fprintf(E->fp,"%s residual (%03d)(%03d) = %.3e from %.3e to %.3e in %5.2f secs \n",
		(convergent ? " * ":"!!!"),count,cycles,residual,r0,acc,CPU_time0()-initial_time);
		fflush(E->fp);
		}
	
    count++;

    free ((void *) D1); 
    free ((void *) Au);
    free ((void *) r);

    if(E->control.verbose)  {
	printf("Total time for %d loops = %g \n",count,CPU_time0()-initial_time);
	for(i=E->mesh.levmax;i>=E->mesh.levmin;i--) {
	    printf("Level %d, total time = %g\n",i,cost_per_level[i]);
	    printf("     total cycles = %d (%g)\n",total_cycles[i],
		   cost_per_level[i]/(1.0e-5+(float)total_cycles[i]));
      }
	printf("projection time = %g\n",E->monitor.cpu_time_on_mg_maps);
	printf("Del sq u solved to accuracy of  %g \n",residual);
    }

    E->control.total_iteration_cycles += count;
    E->control.total_v_solver_calls += 1;
    
    return(convergent);
}

/* =================================
   recursive multigrid function ....
   ================================= */

double multi_grid(E,d1,F,Au,acc,hl)
     struct All_variables *E;
     double *d1;
     double *F;
     double *Au;
     double acc;
     int hl;  /* higher level of two */
{ 
    double residual,AudotAu;
    void project_vector();
    void interp_vector();
    int lev,dlev,ulev,i,j,Vn,Vnmax,ic,cycles;
    double residuaa,alpha,beta;
    void gauss_seidel();
    void element_gauss_seidel();
    void strip_bcs_from_residual();
    void n_assemble_del2_u();

    double vdot(),conj_grad();
    
    FILE *fp;
    char filename[1000];

    const int levmin = E->mesh.levmin; 
    const int levmax = E->mesh.levmax; 

    double time,CPU_time0();
  
    static int been_here = 0;
    static double *res[MAX_LEVELS],*rhs[MAX_LEVELS],*AU[MAX_LEVELS];  
    static double *vel[MAX_LEVELS],*fl[MAX_LEVELS],*del_vel[MAX_LEVELS];
				/* because it's recursive, need a copy at
				    each level */
    if (0==been_here)   { 
	for(i=E->mesh.levmin;i<=E->mesh.levmax;i++) { 
	    vel[i] = (double *)malloc((E->mesh.NEQ[i] + 2)*sizeof(double));
	    res[i] = (double *)malloc((E->mesh.NEQ[i] + 2)*sizeof(double));
	    rhs[i] = (double *)malloc((E->mesh.NEQ[i] + 2)*sizeof(double));
	    fl [i] = (double *)malloc((E->mesh.NEQ[i] + 2)*sizeof(double));
	    del_vel[i]=(double *)malloc((E->mesh.NEQ[i] + 2)*sizeof(double));
	    AU[i] = (double *)malloc((E->mesh.NEQ[i] + 2)*sizeof(double));
	}
    }
    been_here = 1;
 
    Vnmax = E->control.mg_cycle;

    for(j=0;j<E->mesh.NEQ[levmax];j++)
	fl[levmax][j]=F[j];

	/* Project residual onto all the lower levels */
  
    for(lev=levmax;lev>levmin;lev--) { 
	project_vector(E,lev,fl[lev],fl[lev-1],1);
	strip_bcs_from_residual(E,fl[lev-1],lev-1); 
        }

	/* Solve for the lowest level */
	
    cycles = E->control.v_steps_low;
    (void) conj_grad(E,vel[levmin],fl[levmin],AU[levmin],acc*0.001,&cycles,levmin); 


    for(lev=levmin+1;lev<=levmax;lev++) { 
      time=CPU_time0();

                         /* Utilize coarse solution and smooth at this level */
      interp_vector(E,lev-1,vel[lev-1],vel[lev]);
      strip_bcs_from_residual(E,vel[lev],lev);

      for(j=0;j<E->mesh.NEQ[lev];j++)
           rhs[lev][j]=fl[lev][j];

      for(Vn=1;Vn<=Vnmax;Vn++)   {
                                        /*    Downward stoke of the V    */
        for (dlev=lev;dlev>=levmin+1;dlev--)   {

                                      /* Pre-smoothing  */
            cycles=((dlev==levmax)?E->control.v_steps_high:E->control.down_heavy); 
            ic = ((dlev==lev)?1:0);
	    gauss_seidel(E,vel[dlev],rhs[dlev],AU[dlev],0.01,&cycles,dlev,ic);
                                      /* Update residual  */
	    for(i=0;i<E->mesh.NEQ[dlev];i++)
	       res[dlev][i]  = rhs[dlev][i] - AU[dlev][i];

	                              /* Project residual to the lower levels */
	    project_vector(E,dlev,res[dlev],rhs[dlev-1],1);
            strip_bcs_from_residual(E,rhs[dlev-1],dlev-1);
            }

	
                                        /*    Bottom of the V    */
       cycles = E->control.v_steps_low;
       (void) conj_grad(E,vel[levmin],rhs[levmin],AU[levmin],acc*0.001,&cycles,levmin); 

                                        /*    Upward stoke of the V    */
        for (ulev=levmin+1;ulev<=lev;ulev++)   {
            cycles=((ulev==levmax)?E->control.v_steps_high:E->control.up_heavy);

	    interp_vector(E,ulev-1,vel[ulev-1],del_vel[ulev]);
            strip_bcs_from_residual(E,del_vel[ulev],ulev);
            gauss_seidel(E,del_vel[ulev],res[ulev],AU[ulev],0.01,&cycles,ulev,1);

            AudotAu = vdot(E,AU[ulev],AU[ulev],ulev);
            alpha = vdot(E,AU[ulev],res[ulev],ulev)/AudotAu;

	    for(i=0;i<E->mesh.NEQ[ulev];i++)
	       vel[ulev][i] += alpha*del_vel[ulev][i];

            if (ulev ==levmax)
                for(i=0;i<E->mesh.NEQ[ulev];i++)   {
                  res[ulev][i] -= alpha*AU[ulev][i];
                  }


            }
        }
      }
	
    for(j=0;j<E->mesh.NEQ[levmax];j++)   {
	F[j]=res[levmax][j];
	d1[j]=vel[levmax][j];
        }

    residual = sqrt(vselfdot(E,F,levmax)/E->mesh.NEQ[levmax]);

   
    return(residual);
}


/*  ===========================================================
    Conjugate gradient relaxation for the matrix equation Kd = f
    Returns the residual reduction after itn iterations ... 
    ===========================================================  */
 

double conj_grad(E,d0,F,Au,acc,cycles,level)
     struct All_variables *E;
     double *d0;
     double *F;
     double *Au;
     double acc; 
     int *cycles;
     int level;
{ 
    static double *r0,*r1,*r2;
    static double *z0,*z1,*z2;
    static double *p1,*p2;
    static double *Ap;
    static double *BI;
    static double *shuffle;
    static int been_here = 0;

    int count,i,steps;
    double residual;
    double alpha,beta,dotprod,dotr1z1,dotr0z0;

    double CPU_time0(),time;

    void assemble_del2_u();
    void strip_bcs_from_residual();
    double vdot(),vselfdot(),vdot();
    double *dvector();
 
    const int mem_lev=E->mesh.levmax;
    const int high_neq = E->mesh.NEQ[level];
    
    steps = *cycles;

    if (0 == been_here)  /* only used at low level (even if low=high)*/    { 
	r0 = (double *)malloc((1+E->mesh.NEQ[mem_lev])*sizeof(double));
	r1 = (double *)malloc((1+E->mesh.NEQ[mem_lev])*sizeof(double));
	r2 = (double *)malloc((1+E->mesh.NEQ[mem_lev])*sizeof(double));
	z0 = (double *)malloc((1+E->mesh.NEQ[mem_lev])*sizeof(double));
	z1 = (double *)malloc((1+E->mesh.NEQ[mem_lev])*sizeof(double));
	p1 = (double *)malloc((1+E->mesh.NEQ[mem_lev])*sizeof(double));
	p2 = (double *)malloc((1+E->mesh.NEQ[mem_lev])*sizeof(double));
	Ap = (double *)malloc((1+E->mesh.NEQ[mem_lev])*sizeof(double));
	been_here++; 
    }

    for(i=0;i<high_neq;i++)  {
        r1[i] = F[i]; 
        d0[i] = 0.0; 
        }
    
    residual = sqrt(vdot(E,r1,r1,level)/E->mesh.NEQ[level]);

    assert(residual != 0.0  /* initial residual for CG = 0.0 */);
    count = 0;
 
    while (((residual > acc) && (count < steps)) || count == 0)  { 

    for(i=0;i<high_neq;i++)       
	    z1[i] = E->BI[level][i] * r1[i];

	dotr1z1 = vdot(E,r1,z1,level);

	if (0==count)
	    for(i=0;i<high_neq;i++)
	 	  p2[i] = z1[i];                   
	else {
	    assert(dotr0z0 != 0.0 /* in head of conj_grad */);
	    beta = dotr1z1/dotr0z0;
	    for(i=0;i<high_neq;i++)
		  p2[i] = z1[i] + beta * p1[i];
	}

    dotr0z0 = dotr1z1;

	assemble_del2_u(E,p2,Ap,level,1);
    
	dotprod=vdot(E,p2,Ap,level);
	
	if(0.0==dotprod)
	    alpha=1.0e-3;
	else
	    alpha = dotr1z1/dotprod;

	for(i=0;i<high_neq;i++) {
	    d0[i] += alpha * p2[i];
	    r2[i] = r1[i] - alpha * Ap[i];
	}

	residual = sqrt(vdot(E,r2,r2,level)/E->mesh.NEQ[level]);  
	
	shuffle = r0; r0 = r1; r1 = r2; r2 = shuffle;
	shuffle = z0; z0 = z1; z1 = shuffle;
	shuffle = p1; p1 = p2; p2 = shuffle;
		
	count++;
	/* end of while-loop */ 

    }
  
    *cycles=count;
   
    strip_bcs_from_residual(E,d0,level); 
  
   
    return(residual);   }
    

/* ===========================================================================

 */


void jacobi(E,d0,F,Ad,acc,cycles,level,guess)
     struct All_variables *E;
     double *d0;
     double *F,*Ad;
     double acc; 
     int *cycles;
     int level;
     int guess;
{
    static double *r1;
    static double *x1;
    static int been_here=0;
  
    double *dvector();
    void assemble_del2_u();

    int count,steps;
    int i,j,k;

    const int neq=E->mesh.NEQ[level];

    if (0 == been_here) {
	r1 = dvector(0,E->mesh.neq);
	x1 = dvector(0,E->mesh.neq);
	been_here++; 
    }


    steps=*cycles;
    count=0;
   
    if(guess) {
	assemble_del2_u(E,d0,Ad,level,1);
	for(i=0;i<neq;i++)
	    r1[i] = F[i]-Ad[i];
    }
    else
	for(i=0;i<neq;i++)
	    r1[i] = F[i];
   
    while (count <= steps) {
	for(i=0;i<neq;i++) {
	    x1[i] = r1[i] * E->BI[level][i];
	}

	for(i=0;i<neq;i++)
	    r1[i] -= Ad[i];

	count++;
	
    }
    
    for(i=0;i<neq;i++)
	d0[i]=x1[i];

    *cycles=count;
    return;
}

/* ========================================================================================
   An element by element version of the gauss-seidel routine. Initially this is a test 
   platform, we want to know if it handles discontinuities any better than the node/equation
   versions
   =========================================================================================*/

void element_gauss_seidel(E,d0,F,Ad,acc,cycles,level,guess)
    struct All_variables *E;
    double *d0;
    double *F,*Ad;
    double acc; 
    int *cycles;
    int level;
    int guess;
{  
    int count,i,j,k,l,m,ns,nc,d,steps,loc;
    int p1,p2,p3,q1,q2,q3;
    int e,eq,node,node1;
    int element,eqn1,eqn2,eqn3,eqn11,eqn12,eqn13;

    void e_assemble_del2_u();
    void n_assemble_del2_u();
    void strip_bcs_from_residual();
    void get_elt_k();
    
    double U1[24],AD1[24],F1[24];
    double w1,w2,w3;
    double w11,w12,w13;
    double w[24];
    static double *Ad0,*dd,*elt_k;
    static int *vis,been_here=0;

    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[E->mesh.nsd];
    const int neq=E->mesh.NEQ[level];
    const int nel=E->mesh.NEL[level];
    const int nno=E->mesh.NNO[level];

   
    steps=*cycles;
   
    if(0==been_here) {
	dd = (double *)malloc((neq+1)*sizeof(double));
	vis = (int *)malloc((nno+1)*sizeof(int));
	elt_k=(double *)malloc((24*24)*sizeof(double));
	been_here++;
    }
    
    if(guess){
	e_assemble_del2_u(E,d0,Ad,level,1);
    }
    else {
	for(i=0;i<neq;i++) 
	     Ad[i]=d0[i]=0.0;
    }
	
    count=0;
    while (count <= steps) {
	for(i=1;i<=nno;i++)
	    vis[i]=0;
	
	for(e=1;e<=nel;e++) {

	    elt_k = E->elt_k[level][e].k;
	    
	    for(i=1;i<=ends;i++)  {
		node=E->IEN[level][e].node[i];
		p1=(i-1)*dims; 
		w[p1] = w[p1+1] = 1.0;
		if(dims==3)
		    w[p1+2] = 1.0;
		if(E->NODE[level][node] & VBX)
		    w[p1] = 0.0;
		if(E->NODE[level][node] & VBZ)
		    w[p1+1] = 0.0;
		if((3==dims) && (E->NODE[level][node] & VBY))
		    w[p1+2] = 0.0;
		
	    }


	    for(i=1;i<=ends;i++){
		node=E->IEN[level][e].node[i];
		if(!vis[node])
		    continue;
		
		eqn1=E->ID[level][node].doff[1];
		eqn2=E->ID[level][node].doff[2];
		if(3==dims)
		    eqn3=E->ID[level][node].doff[3];
		p1=(i-1)*dims*n;
		p2=p1+n; 
		p3=p2+n;

	
		/* update Au */
		for(j=1;j<=ends;j++) {
		    node1=E->IEN[level][e].node[j];
		   
		    if(3==dims){
			eqn11=E->ID[level][node1].doff[1];
			eqn12=E->ID[level][node1].doff[2];
		   	eqn13=E->ID[level][node1].doff[3];
			q1=(j-1)*3;
		
			Ad[eqn11] += w[q1]*(elt_k[p1+q1] * dd[eqn1] + elt_k[p2+q1] * dd[eqn2] + elt_k[p3+q1] * dd[eqn3]);
		   	Ad[eqn12] += w[q1+1]*(elt_k[p1+q1+1] * dd[eqn1] + elt_k[p2+q1+1] * dd[eqn2] + elt_k[p3+q1+1] * dd[eqn3]);
			Ad[eqn13] += w[q1+2]*(elt_k[p1+q1+2] * dd[eqn1] + elt_k[p2+q1+2] * dd[eqn2] + elt_k[p3+q1+2] * dd[eqn3]);
		    }  
		    
		    else {
			eqn11=E->ID[level][node1].doff[1];
			eqn12=E->ID[level][node1].doff[2];
			q1=(j-1)*2;
		
			Ad[eqn11] += w[q1] * (elt_k[p1+q1] * dd[eqn1] + elt_k[p2+q1] * dd[eqn2]);
			Ad[eqn12] += w[q1+1] * (elt_k[p1+q1+1] * dd[eqn1] + elt_k[p2+q1+1] * dd[eqn2]);
		    }		    
		}
	    } 

 
	    for(i=1;i<=ends;i++){
		node=E->IEN[level][e].node[i];
		if(vis[node])
		    continue;

		eqn1=E->ID[level][node].doff[1];
		eqn2=E->ID[level][node].doff[2];
		if(3==dims) {
		    eqn3=E->ID[level][node].doff[3];
		}
		p1=(i-1)*dims*n;
		p2=p1+n; 
		p3=p2+n;
		
		/* update dd, d0 */
		d0[eqn1] += (dd[eqn1] = w[(i-1)*dims]*(F[eqn1]-Ad[eqn1])*E->BI[level][eqn1]); 
		d0[eqn2] += (dd[eqn2] = w[(i-1)*dims+1]*(F[eqn2]-Ad[eqn2])*E->BI[level][eqn2]); 
		if(3==dims){
		    d0[eqn3] += (dd[eqn3] = w[(i-1)*dims+2]*(F[eqn3]-Ad[eqn3])*E->BI[level][eqn3]); 
		}
		
		vis[node]=1;
		
		/* update Au */
		for(j=1;j<=ends;j++) {
		   node1=E->IEN[level][e].node[j];
		    
		   if(3==dims){
		       eqn11=E->ID[level][node1].doff[1];
		       eqn12=E->ID[level][node1].doff[2];
		       eqn13=E->ID[level][node1].doff[3];
		       q1=(j-1)*3;
		       q2=q1+1;
		       q3=q1+2;
		       
		       Ad[eqn11] += w[q1]*(elt_k[p1+q1] * dd[eqn1] + elt_k[p2+q1] * dd[eqn2] + elt_k[p3+q1] * dd[eqn3]);
		       Ad[eqn12] += w[q2]*(elt_k[p1+q2] * dd[eqn1] + elt_k[p2+q2] * dd[eqn2] + elt_k[p3+q2] * dd[eqn3]);
		       Ad[eqn13] += w[q3]*(elt_k[p1+q3] * dd[eqn1] + elt_k[p2+q3] * dd[eqn2] + elt_k[p3+q3] * dd[eqn3]);
		   }  
		   else {
		       eqn11=E->ID[level][node1].doff[1];
		       eqn12=E->ID[level][node1].doff[2];
		       q1=(j-1)*2;
		       q2=q1+1;
		       Ad[eqn11] += w[q1] * (elt_k[p1+q1] * dd[eqn1] + elt_k[p2+q1] * dd[eqn2]);
		       Ad[eqn12] += w[q2] * (elt_k[p1+q2] * dd[eqn1] + elt_k[p2+q2] * dd[eqn2]);
		   }		   
		   
		}
	    }
	
	}
	/* completed cycle */

	    count++;
	
    }
    
    return;
}


/* ============================================================================
   Multigrid Gauss-Seidel relaxation scheme which requires the storage of local
   information, otherwise some other method is required. NOTE this is a bit worse
   than real gauss-seidel because it relaxes all the equations for a node at one
   time (Jacobi at a node). It does the job though.
   ============================================================================ */

void gauss_seidel(E,d0,F,Ad,acc,cycles,level,guess)  
     struct All_variables *E;
     double *d0;
     double *F,*Ad;
     double acc; 
     int *cycles;
     int level;
     int guess;
{  
   
    int count,i,j,k,l,m,ns,steps;
    int *C;
    int eqn1,eqn2,eqn3,max_eqn;
    
    void n_assemble_del2_u();
    double vdot(),vselfdot();
    double U1,U2,U3;
   
    higher_precision *B1,*B2,*B3;
  
    
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[E->mesh.nsd];
    const int neq=E->mesh.NEQ[level];
    const int num_nodes=E->mesh.NNO[level];
    const int nox=E->mesh.NOX[level];
    const int noz=E->mesh.NOY[level];
    const int noy=E->mesh.NOZ[level];

    steps=*cycles;
  	       
    if(guess) {
	d0[neq]=0.0;
	n_assemble_del2_u(E,d0,Ad,level,1);
    }
    else
	for(i=0;i<neq;i++) {
	    d0[i]=Ad[i]=0.0;
	}

    count = 0;

 /* fprintf(E->fp,"in g_s level=%d\n",level);
  for(i=0;i<neq;i++)
    fprintf(E->fp,"%d %g %g\n",i,F[i],E->BI[level][i]);
*/ 
    while (count < steps) {
	for(m=1;m<=noy;m++)
	    for(l=1;l<=noz;l++)  
		for(k=1;k<=nox;k++)
		{
		    i=l+(k-1)*noz+(m-1)*noz*nox;

	    if(E->NODE[level][i] & OFFSIDE) 
		continue;

  	    max_eqn = E->Node_eqn[level][i+1]-E->Node_eqn[level][i];
            C=E->Node_map[level]+E->Node_eqn[level][i];
	    B1=E->Eqn_k[level]+E->Node_k_id[level][i];
	    B2=E->Eqn_k[level]+E->Node_k_id[level][i]+max_eqn;
	   		
	    if(2==dims){
		eqn1=E->ID[level][i].doff[1];
		eqn2=E->ID[level][i].doff[2];
		U1 = (F[eqn1] - Ad[eqn1])*E->BI[level][eqn1]; 
		U2 = (F[eqn2] - Ad[eqn2])*E->BI[level][eqn2];
		d0[eqn1] += U1;
		d0[eqn2] += U2;
/*    fprintf(E->fp,"%d %d %d %g %g %g %g\n",max_eqn,eqn1,eqn2,U1,U2,Ad[eqn1],Ad[eqn2]);
*/		    
		for(j=0;j<max_eqn;j++) {
		    Ad[C[j]]  += B1[j] * U1 +  B2[j] * U2; 
/*    fprintf(E->fp,"%d %d %g %g \n",j,C[j],B1[j],B2[j]);  */
		} 
	    }
	    else {
		eqn1=E->ID[level][i].doff[1];
		eqn2=E->ID[level][i].doff[2];
		eqn3=E->ID[level][i].doff[3];
		B3=E->Eqn_k[level]+E->Node_k_id[level][i]+2*max_eqn;
		U1 = (F[eqn1] - Ad[eqn1])*E->BI[level][eqn1]; 
		U2 = (F[eqn2] - Ad[eqn2])*E->BI[level][eqn2];
		U3 = (F[eqn3] - Ad[eqn3])*E->BI[level][eqn3];
		d0[eqn1] += U1;
		d0[eqn2] += U2;
		d0[eqn3] += U3;
		   
		for(j=0;j<max_eqn;j++) {
		    Ad[C[j]]  += B1[j]*U1 +  B2[j]*U2 +  B3[j]*U3; 
		} 
	    }
		
		
	}
	    
	count++;
    }
    
       

    *cycles=count;
    return;
    
}
 
/* ==========================================================================
   This is the full gauss-seidel method. Slightly more accuracy is gained by
   relaxing over every equation in turn but the access to memory is probably
   better in the routine above
   ========================================================================== */

void gauss_seidel1(E,d0,F,Ad,acc,cycles,level,guess)  
     struct All_variables *E;
     double *d0;
     double *F,*Ad;
     double acc; 
     int *cycles;
     int level;
     int guess;
{  
    static double *r1;
    static int been_here = 0;
  
    int count,i,j,ns,steps,max_eqn;
    int *C;
    int eqn1,eqn2,eqn3;
    
    float CPU_time();

    void assemble_del2_u();
    void n_assemble_del2_u();
    double vdot(),vselfdot();
    double A,U1,U2,U3;
    double residual,r0,R0;
   
    higher_precision *B1,*B2,*B3;
  
    
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[E->mesh.nsd];
    const int neq=E->mesh.NEQ[level];
    const int num_nodes=E->mesh.NNO[level];
 
   steps=*cycles;
   
    if(0==been_here){
	r1=(double *)malloc((E->mesh.neq+2)*sizeof(double));
	been_here++;
    }

 
	if(guess) {
	    d0[neq]=0.0;
	    n_assemble_del2_u(E,d0,Ad,level,1);
	    for(i=0;i<neq;i++) {
		r1[i]=F[i]-Ad[i];
	    }
	}
	else
	    for(i=0;i<neq;i++) {
		r1[i]=F[i];
		d0[i]=Ad[i]=0.0;
	    }

	/* initial residual */
	R0=1.5*(residual=sqrt(vselfdot(E,r1,level)/neq));
	r0=1.0/residual;

	count = 1;
	while (count <= steps && residual*r0 > acc && residual <= R0) {
	    for(i=1;i<=num_nodes;i++) {
		if(E->NODE[level][i] & OFFSIDE) 
		    continue;
		
		eqn1=E->ID[level][i].doff[1];
		eqn2=E->ID[level][i].doff[2];

  	        max_eqn = E->Node_eqn[level][i+1]-E->Node_eqn[level][i];
                C=E->Node_map[level]+E->Node_eqn[level][i];
	        B1=E->Eqn_k[level]+E->Node_k_id[level][i];
	        B2=E->Eqn_k[level]+E->Node_k_id[level][i]+max_eqn;

		assert(eqn1 >= 0 && eqn1 < neq);
		assert(eqn2 >= 0 && eqn2 < neq);
		 
		if(3==dims) {
		    eqn3=E->ID[level][i].doff[3];
	            B3=E->Eqn_k[level]+E->Node_k_id[level][i]+2*max_eqn;
		    assert(eqn3 >= 0 && eqn3 < neq);
		}

		/* dirn 1 */

		U1 = (F[eqn1] - Ad[eqn1]) * E->BI[level][eqn1];
		d0[eqn1] += U1;
		for(j=0;j<max_eqn;j++) {
		    Ad[C[j]-1]  += B1[j] * U1; 
		}
  
		/* dirn 2 */
		U2 = (F[eqn2] - Ad[eqn2]) * E->BI[level][eqn2];
		d0[eqn2] +=  U2;
		for(j=0;j<max_eqn;j++) {
		    Ad[C[j]-1]  += B2[j] * U2; 
		}

		/* dirn 3 */
		if(3==dims) {
		    U3 = (F[eqn3] - Ad[eqn3]) * E->BI[level][eqn3];
		    d0[eqn3] += U3; 
		    for(j=0;j<max_eqn;j++) {
			Ad[C[j]-1]  += B3[j] * U3; 
		    }
		}
	    }
		

	    for(i=0;i<neq;i++) {
		r1[i]=F[i]-Ad[i];
	    }
	    residual=sqrt(vselfdot(E,r1,level)/neq);

	    count++;
	} 

    *cycles=count-1;
    return;
    
}
 
void print_elt_k(E,a)
     struct All_variables *E;
     double a[24*24];

{ int l,ll,n;
 
  printf("elt k is ...\n");

 
  n = loc_mat_size[E->mesh.nsd];
  
  for(l=0;l<n;l++)
    { fprintf(stderr,"\n");fflush(stderr);
      for(ll=0;ll<n;ll++)
	{ fprintf(stderr,"%s%.3e ",a[ll*n+l] >= 0.0 ? "+" : "",a[ll*n+l]);
	  fflush(stderr);
	}
    }
  fprintf(stderr,"\n"); fflush(stderr);

  return; }


double cofactor(A,i,j,n)
     double A[4][4];
     int i,j,n;

{ int k,l,p,q;
  double determinant();
  static int been_here = 0;
  double B[4][4]; /* because of recursive behaviour of det/cofac, need to use
			       new copy of B at each 'n' level of this routine */
 
  if (n>3) printf("Error, no cofactors for matrix more than 3x3\n");

  p=q=1;

  for(k=1;k<=n;k++)    {
     if(k==i) continue;
     for(l=1;l<=n;l++)      {
	   if (l==j) continue;
           B[p][q]=A[k][l];
	   q++ ;
	   }
     q=1;p++;  
     }

       
  return(epsilon[i][j]*determinant(B,n-1));

 
}


/* Fast (conditional) determinant for 3x3 or 2x2 ... otherwise calls general routine */

double determinant(A,n)
     double A[4][4];
     int n;

{ double gen_determinant();

  switch (n)
    { case 1: 
	return(A[1][1]);
	break;
      case 2:
	return(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	break;
      case 3:
	return(A[1][1]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])-
	       A[1][2]*(A[2][1]*A[3][3]-A[2][3]*A[3][1])+
	       A[1][3]*(A[2][1]*A[3][2]-A[2][2]*A[3][1]));
	break;
      default:
	return(1);
/*	return(gen_determinant(A,n)); */
      }
}
			 

/* recursive function to determine matrix determinant */

double gen_determinant(A,n)
     double **A;
     int n;

{ double det;
  double cofactor();

  int i;
 
 
  if(n==1) return(A[1][1]); /* need a way to break the recursion */
  
  det=0.0; 
  for(i=1;i<=n;i++)
    det += A[1][i]*cofactor(A,1,i,n);
     
  return(det);
}


 float area_of_4node(x1,y1,x2,y2,x3,y3,x4,y4)
 float x1,y1,x2,y2,x3,y3,x4,y4;

 {
 float area;

 area = fabs(0.5*(x1*(y2-y4)+x2*(y4-y1)+x4*(y1-y2)))
      + fabs(0.5*(x2*(y3-y4)+x3*(y4-y2)+x4*(y2-y3)));

 return area;
 }


float Tmax(E,T)
  struct All_variables *E;
  float *T;
{
  float temp;
  int i;

  temp = -10.0;
  for(i=1;i<=E->mesh.nno;i++)
    temp = max(T[i],temp);

  return (temp);
  }


float  vnorm_nonnewt(E,dU,U,lev)
  struct All_variables *E;
  float *dU,*U;
  int lev;
{
 float temp1,temp2,dtemp,temp;
 int a,e,i;
 const int dims = E->mesh.nsd;
 const int ends = enodes[dims];
 const int nel=E->mesh.nel;

 dtemp=0.0;
 temp=0.0;
 for (e=1;e<=nel;e++)
   if (E->mat[e]==1)
     for (i=1;i<=dims;i++)
       for (a=1;a<=ends;a++) {
         dtemp += dU[ E->LMD[lev][e].node[a].doff[i] ]*
                  dU[ E->LMD[lev][e].node[a].doff[i] ];
         temp +=  U[ E->LMD[lev][e].node[a].doff[i] ]*
                  U[ E->LMD[lev][e].node[a].doff[i] ];
         }

  temp1 = sqrt(dtemp/temp);

  return (temp1);
  }

/* ===============================================
   strips horizontal average from nodal field X.
   Assumes orthogonal mesh, otherwise, horizontals
   aren't & another method is required.
   =============================================== */

void remove_horiz_ave(E,X,H,store_or_not)
     struct All_variables *E;
     float *X, *H;
     int store_or_not;

{
    int i,j,k,n,ln,nox,noz,noy;
    void return_horiz_ave();

    const int dims = E->mesh.nsd;

    noy = E->mesh.noy;
    noz = E->mesh.noz;
    nox = E->mesh.nox;

    return_horiz_ave(E,X,H);
    for(i=1;i<=noz;i++)
      for(k=1;k<=noy;k++)  
        for(j=1;j<=E->mesh.nox;j++)     {
          n = i+(j-1)*noz+(k-1)*nox*noz;
          X[n] -= H[i];
          }

   return;
  }

void return_horiz_ave(E,X,H)
     struct All_variables *E;
     float *X, *H;
{
  const int dims = E->mesh.nsd;
  int i,j,k,d,nint,noz,nox,noy,el,elz,elx,ely,j1,j2,i1,i2,k1,k2;
  int lnode[5], sizeofH, noz2,iroot;
  float *Have,*temp;

  sizeofH = (2*E->mesh.noz+2)*sizeof(float);

  Have = (float *)malloc(sizeofH);
  temp = (float *)malloc(sizeofH);
  noz = E->mesh.noz;
  noy = E->mesh.noy;
  elz = E->mesh.elz;
  ely = E->mesh.ely;
  elx = E->mesh.elx;
  noz2 = 2*noz;

  for (i=1;i<=elz;i++)  {
    temp[i] = temp[i+noz] = 0.0;
    temp[i+1] = temp[i+1+noz] = 0.0;
    for (j=1;j<=elx;j++)     {
            el = i + (j-1)*elz;

            lnode[1] = E->ien[el].node[1];
            lnode[2] = E->ien[el].node[4];
            temp[i] += 0.5*(X[lnode[1]]+X[lnode[2]])*sin((E->X[1][lnode[2]]+E->X[1][lnode[1]])*0.5)*(E->X[1][lnode[2]]-E->X[1][lnode[1]]);
            temp[i+noz] += sin((E->X[1][lnode[2]]+E->X[1][lnode[1]])*0.5)*(E->X[1][lnode[2]]-E->X[1][lnode[1]]);

        if (i==elz)  {
              lnode[1] = E->ien[el].node[2];
              lnode[2] = E->ien[el].node[3];

              temp[i+1] += 0.5*(X[lnode[1]]+X[lnode[2]])*sin((E->X[1][lnode[2]]+E->X[1][lnode[1]])*0.5)*(E->X[1][lnode[2]]-E->X[1][lnode[1]]);
              temp[i+1+noz] += sin((E->X[1][lnode[2]]+E->X[1][lnode[1]])*0.5)*(E->X[1][lnode[2]]-E->X[1][lnode[1]]);

          }   /* end of if i==elz    */
      }        /* end of j */

     }        /* Done for i */

    for (i=1;i<=noz2;i++)
      Have[i] = temp[i];

  for (i=1;i<=noz;i++) {
    if(Have[i+noz] != 0.0)
       H[i] = Have[i]/Have[i+noz];
    }

  free ((void *) Have);
  free ((void *) temp);
  return;
  }


void return_horiz_ave1(E,X,H)
     struct All_variables *E;
     float *X, *H;
{
  const int dims = E->mesh.nsd;
  int i,j,k,d,nint,noz,nox,noy,el,elz,elx,ely,j1,j2,i1,i2,k1,k2;
  int lnode[5], sizeofH, noz2,iroot;
  float *Have,*temp;
  struct Shape_function1 M;
  struct Shape_function1_dA dGamma;
  void get_global_1d_shape_fn();

  sizeofH = (2*E->mesh.noz+2)*sizeof(float);

  Have = (float *)malloc(sizeofH);
  temp = (float *)malloc(sizeofH);

  noz = E->mesh.noz;
  noy = E->mesh.noy;
  elz = E->mesh.elz;
  ely = E->mesh.ely;
  elx = E->mesh.elx;
  noz2 = 2*noz;

  for (i=1;i<=elz;i++)  {
    temp[i] = temp[i+noz] = 0.0;
    temp[i+1] = temp[i+1+noz] = 0.0;
    for (k=1;k<=ely;k++)    
      for (j=1;j<=elx;j++)     {
            el = i + (j-1)*elz + (k-1)*elx*elz;
            get_global_1d_shape_fn(E,el,&M,&dGamma,0);

            lnode[1] = E->ien[el].node[1];
            lnode[2] = E->ien[el].node[4];
            lnode[4] = E->ien[el].node[5];
            lnode[3] = E->ien[el].node[8];

        for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
          for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
              temp[i] += X[lnode[d]] * E->M.vpt[GMVINDEX(d,nint)]
                         * dGamma.vpt[GMVGAMMA(1,nint)];
              temp[i+noz] += E->M.vpt[GMVINDEX(d,nint)]
                         * dGamma.vpt[GMVGAMMA(1,nint)];
              }

        if (i==elz)  {
              lnode[1] = E->ien[el].node[2];
              lnode[2] = E->ien[el].node[3];
              lnode[4] = E->ien[el].node[6];
              lnode[3] = E->ien[el].node[7];

          for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
            for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
              temp[i+1] += X[lnode[d]] * E->M.vpt[GMVINDEX(d,nint)]
                         * dGamma.vpt[GMVGAMMA(1,nint)];
              temp[i+1+noz] += E->M.vpt[GMVINDEX(d,nint)]
                         * dGamma.vpt[GMVGAMMA(1,nint)];
              }
          }   /* end of if i==elz    */
      }        /* end of k */

     }        /* Done for i */

    for (i=1;i<=noz2;i++)
      Have[i] = temp[i];

  for (i=1;i<=noz;i++) {
    if(Have[i+noz] != 0.0)
       H[i] = Have[i]/Have[i+noz];
    }

  free ((void *) Have);
  free ((void *) temp);
  return;
  }

float return_bulk_value(E,Z,average)
     struct All_variables *E;
     float *Z;
     int average;

{ 
    void get_global_shape_fn();
    void float_global_operation();

    double xk[3][5];
    int n,i,j,k,el;
    float volume,integral,volume1,integral1;
  
    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;

    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];

    volume=0.0;
    integral=0.0;
    for (el=1;el<=E->mesh.nel;el++)  {

          get_global_shape_fn(E,el,&GN,&GNx,&dOmega,xk,0,E->mesh.levmax);

          for(j=1;j<=vpts;j++)
            for(i=1;i<=ends;i++) {
                n = E->ien[el].node[i];
                volume += E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
                integral += Z[n] * E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
            }

      }


    if(average && volume != 0.0)
           integral /= volume;

    return((float)integral);
   }

double kineticE(E,A,lev)
   struct All_variables *E;
   double *A;
   int lev;

{
  int i,neq;
  double temp;

  temp = 0.0;

  neq=E->mesh.NEQ[lev];
  for (i=0;i<neq;i++)
      temp += A[i]*A[i];

  temp = temp/neq;
  return (temp);
}

/****************************************************************
 * gaussj																												*
 *																															*
 * Linear equation solution by Gauss-Jordan elimination, 				*
 * a[1..n][1..n] is the input matrix.  b[1..n] is input					*
 * containing the right-hand side vecotr.  On output, a is 			*
 * replaced by its matrix inverse, and b is replaced by the 		*
 * corresponding set of solution vectors.												*
 *																															*
 * Note: Adapted from Numerical Recipes													*
 *																															*
 * Parameters: a	n by n matrix																	*
 *						 n	length of vectors, size of matrix							*
 *						 b  rhs n-vector																	*
 *																															*
 ****************************************************************/
void gaussj(E,a,n,b)
	struct All_variables *E;
	complex double **a;
	int n;
	complex double *b;
{
	int *indxc,*indxr,*ipiv;		/* For bookkeeping on the pivoting */
	int i,icol,irow,j,k,l,ll;		/* counters */
	double big;
	complex double pivinv,dum,temp;

	indxc = (int *)malloc((n+1)*sizeof(int));
	indxr = (int *)malloc((n+1)*sizeof(int));
	ipiv = (int *)malloc((n+1)*sizeof(int));
	
	for (j=1;j<=n;j++)
		ipiv[j] = 0;

	for (i=1;i<=n;i++) {				/* Main loop over the columns to be reduced */
		big=0.0;
		for (j=1;j<=n;j++) 			/* Outer loop of the search for a pivot element */
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (cabs(a[j][k]) >= big) {
							big=cabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} 
					else if (ipiv[k] > 1)
						fprintf(E->fp,"gaussj: Singular Matrix-1\n");
				}
		++ipiv[icol];

		/* 
		 *	We now have the pivot element, so we interchange rows if needed, 
		 *	to put the pivot element on the diagonal.  The columns are not 
		 *	physically interchanged, only relabeled: indxc[i], the column of
		 *	the ith pivot element, is the ith column that is reduced, while
		 *	indxr[i] is the row in which that pivot element was originally 
		 *	located.  If indxr[i] != indxc[i], there is an implied column 
		 *	interchange.  With this form of bookkeeping, the solution b's will
		 *	end up in the correct order, and the inverse matrix will be
		 *	scrambled by columns.
		 */

		if (irow != icol) {
			for(l=1;l<=n;l++)
				SWAP(a[irow][l],a[icol][l]);
			SWAP(b[irow],b[icol]);	/* b is a vector */
		}

		/* We are now ready to divide the pivot row by the pivot element, 
			 located at irow and icol. */
		indxr[i] = irow;				
		indxc[i] = icol;				
		if (a[icol][icol] == 0.0)
			fprintf(E->fp,"gaussj: Singular Matrix-2\n");
		pivinv = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for(l=1;l<=n;l++)
			a[icol][l] *= pivinv;
		b[icol] *= pivinv;

		/* Next we reduce the rows, except for the pivot on, of course. */
		for(ll=1; ll<=n; ll++)
			if(ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for(l=1;l<=n;l++)
					a[ll][l] -= a[icol][l]*dum;
				b[ll] -= b[icol]*dum;
			}
	}

	/* 
	 * This is the end of the main loop over columns of the reduction.
	 * It only remains to unscramble the solution in view of the column 
	 * interchanges.  We do this by interchanging pairs of columns in
	 * the reverse order that the permutation was built up. 
	 */
	for(l=n;l>=1;l--) {
		if(indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}

	/* And we are done */
  free((void *) ipiv); 
  free((void *) indxr); 
  free((void *) indxc); 
}

