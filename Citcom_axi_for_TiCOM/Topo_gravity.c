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

#define c_re(a) a.real
#define c_im(a) a.imag
  typedef struct compl {  double real;
			  double imag;    } COMPLEX;

extern float Zj0[1000],Zj1[1000];


/* ===================================================================
   Consistent boundary flux method for stress ... Zhong,Gurnis,Hulbert 

   Solve for the stress as the code defined it internally, rather than
   what was intended to be solved. This is more appropriate.

   Note also that the routine is dependent on the method 
   used to solve for the velocity in the first place.
   ===================================================================  */



void get_CBF_topo(E,H,HB)       /* call this only for top and bottom processors*/
    struct All_variables *E;
    float *H,*HB;
   
{
    void get_elt_k();
    void get_elt_g();
    void get_elt_f();
    void matrix_transform_g();
    void get_global_1d_shape_fn();

    int el,elb,els,node,nodeb,nodes,i,j,k,l,m,n,count;
    int nodel,nodem,nodesl,nodesm,lnsf,nel2;
    
    struct Shape_function1 GM,GMb;
    struct Shape_function1_dA dGammax,dGammabx;
 
    float *eltTU,*eltTL,*SU,*SL,*RU,*RL;

    double eltk[24*24],eltf[24];
    double eltkb[24*24],eltfb[24];
    double res[24],resb[24],eu[24],eub[24];
    higher_precision eltg[24][1],eltgb[24][1];

    const int dims=E->mesh.nsd;
    const int Tsize=5;   /* maximum values, applicable to 3d, harmless for 2d */ 
    const int Ssize=4;
    const int ends=enodes[dims];
    const int noz=E->mesh.noz;
    const int noy=E->mesh.noy;
    const int nno=E->mesh.nno;
    const int onedv=onedvpoints[dims];
    const int snode1=1,snode2=4,snode3=5,snode4=8;
    const int elz = E->mesh.elz;
    const int ely = E->mesh.ely;
    const int lev=E->mesh.levmax;
 
    lnsf=E->mesh.nsf;
 
    eltTU = (float *)malloc((1+Tsize)*sizeof(float)); 
    eltTL = (float *)malloc((1+Tsize)*sizeof(float));
    SU = (float *)malloc((1+lnsf)*sizeof(float));
    SL = (float *)malloc((1+lnsf)*sizeof(float));
    RU = (float *)malloc((1+lnsf)*sizeof(float));
    RL = (float *)malloc((1+lnsf)*sizeof(float));

    for(i=0;i<=lnsf;i++)
      RU[i] = RL[i] = SU[i] = SL[i] = 0.0;

    /* calculate the element residuals */

    for(els=1;els<=E->mesh.snel;els++) {
	    el = E->surf_element[els];
	    elb = el + elz-1;

	    for(m=0;m<ends;m++) {  /* for bottom elements */
          nodeb= E->ien[elb].node[m+1];
          eub[m*dims  ] = E->V[1][nodeb];
          eub[m*dims+1] = E->V[2][nodeb];
          if(3==dims) 
            eub[m*dims+2] = E->V[3][nodeb]; 
          }

	      for(m=0;m<ends;m++) {  
          node = E->ien[el].node[m+1];
          eu [m*dims  ] = E->V[1][node];
          eu [m*dims+1] = E->V[2][node];
          if(3==dims)
            eu [m*dims+2] = E->V[3][node];
          }

	    get_elt_k(E,el,eltk,lev,0);
	    get_elt_k(E,elb,eltkb,lev,0);
	    get_elt_f(E,el,eltf,0,0);
	    get_elt_f(E,elb,eltfb,0,0);
        get_elt_g(E,el,eltg,lev);
	    get_elt_g(E,elb,eltgb,lev);
	   
	    for(m=0;m<dims*ends;m++) {
          res[m]  = eltf[m]  - E->elt_del[lev][el].g[m][0]  * E->P[el];
          resb[m] = eltfb[m] - E->elt_del[lev][elb].g[m][0]* E->P[elb];
	      }
	   
	    for(m=0;m<dims*ends;m++)
		for(l=0;l<dims*ends;l++) {
		    res[m]  -= eltk[ends*dims*m+l]  * eu[l];
		    resb[m] -= eltkb[ends*dims*m+l] * eub[l];
		   }
	   
	    /* Put relevant (vertical & surface) parts of element residual into surface residual */
		
	    for(m=1;m<=ends;m++) {     /* for bottom elements */
		switch (m) {
		case 2:
		    RL[E->sien[els].node[1]] += resb[(m-1)*dims+1];  
		    break;
		case 3:
		    RL[E->sien[els].node[2]] += resb[(m-1)*dims+1];  
		    break;
		case 7:
		    RL[E->sien[els].node[3]] += resb[(m-1)*dims+1];  
		    break;
		case 6:
		    RL[E->sien[els].node[4]] += resb[(m-1)*dims+1];  
		    break;
		    }
		}


	    for(m=1;m<=ends;m++) {
		switch (m) {
		case 1:
		    RU[E->sien[els].node[1]] += res[(m-1)*dims+1];  
		    break;
		case 4:
		    RU[E->sien[els].node[2]] += res[(m-1)*dims+1];  
		    break;
		case 8:
		    RU[E->sien[els].node[3]] += res[(m-1)*dims+1];  
		    break;
		case 5:
		    RU[E->sien[els].node[4]] += res[(m-1)*dims+1];  
		    break;
		    }
		}
	}
    
    /* calculate the LHS */
 
    for(els=1;els<=E->mesh.snel;els++) {
	    el = E->surf_element[els];
	    elb = el + elz-1;

	    get_global_1d_shape_fn(E,el,&GM,&dGammax,1);
	    get_global_1d_shape_fn(E,elb,&GMb,&dGammabx,1);
   
	    for(m=1;m<=onedv;m++)        {
	      eltTU[m-1] = 0.0;
	      eltTL[m-1] = 0.0; 
	      for(n=1;n<=onedv;n++)          {
		     eltTU[m-1] += 
		         dGammax.vpt[GMVGAMMA(1,n)] * l_1d[n].weight[dims-1]
			 * E->L.vpt[GMVINDEX(m,n)] * E->L.vpt[GMVINDEX(m,n)];
		     eltTL[m-1] += 
			 dGammabx.vpt[GMVGAMMA(1+dims,n)]*l_1d[n].weight[dims-1]
			 * E->L.vpt[GMVINDEX(m,n)] * E->L.vpt[GMVINDEX(m,n)];
		     }
		  }

        for (m=1;m<=onedv;m++)     /* for bottom */
	      SL[E->sien[els].node[m]] += eltTL[m-1];

          for (m=1;m<=onedv;m++) 
	           SU[E->sien[els].node[m]] += eltTU[m-1];
	    }


      for(i=1;i<=E->mesh.nsf;i++)
        H[i] = -RU[i]/SU[i];
        
      for(i=1;i<=E->mesh.nsf;i++)
        HB[i] = -RL[i]/SL[i];

    free((void *)eltTU);
    free((void *)eltTL);
    free((void *)SU);
    free((void *)SL);
    free((void *)RU);
    free((void *)RL);
    return;

}

void get_STD_topo(E,tpg,tpgb,ii)
    struct All_variables *E;
    float *tpg;
    float *tpgb;
    int ii;
{
    int i,j,k,snode,node;
    float *Szz;
    void get_Szz();

    Szz=(float *) malloc((1+E->mesh.nno)*sizeof(float));
    
    get_Szz(E,E->P,Szz,ii);

    for(snode=1;snode<=E->mesh.nsf;snode++)   {
           node = E->surf_node[snode];
	   tpg[snode]  = -2*Szz[node] + Szz[node-1];
	   tpgb[snode] = 2*Szz[node-E->mesh.noz+1]-Szz[node-E->mesh.noz+2]; 
	   }

    free((void *)Szz);
    return;
}

void get_Szz(E,P,SZZ,file_number)
     struct All_variables *E;
     double *P;
     int file_number;
     float *SZZ;

{
    void get_surf_stress();
    void get_global_shape_fn();
    int i,j,k,e,node, nel2;
    
    float *SXX,*SYY,*SXY,*SXZ,*SZY,VZ[9],VY[9],VX[9],Szz,Sxx,Syy,Sxy,Sxz,Szy;
    float EXX[9],EZZ[9],EXZ[9],EYY[9];
    double el_volume,visc[9],Visc,a,b,xk[3][5];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;

/*
    SXX = (float *)malloc((E->mesh.nno+1)*sizeof(float));
    SYY = (float *)malloc((E->mesh.nno+1)*sizeof(float));
    SXY = (float *)malloc((E->mesh.nno+1)*sizeof(float));
    SXZ = (float *)malloc((E->mesh.nno+1)*sizeof(float));
    SZY = (float *)malloc((E->mesh.nno+1)*sizeof(float));
*/

    for(i=1;i<=E->mesh.nno;i++) {
      SZZ[i] = 0.0;
      }

    for(e=1;e<=E->mesh.nel;e++)  {
      Szz = 0.0;
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,xk,0,E->mesh.levmax);
	
      for(i=1;i<=vpts;i++)   {
	  visc[i] =  E->EVI[E->mesh.levmax][(e-1)*vpts+i] * dOmega.vpt[i];
          EZZ[i] = 0.0;
          }
     
      for(j=1;j<=ends;j++) {
          VZ[j] = E->V[2][E->ien[e].node[j]];
	  }

      for(i=1;i<=vpts;i++)  {
        for(j=1;j<=ends;j++)  {
          EZZ[i] += VZ[j]* GNx.vpt[GNVXINDEX(1,j,i)];
	  }
        Szz += 2.0 * visc[i] * EZZ[i];
	}

      Szz /= E->eco[e].area;

      Szz -= P[e];  /* add the pressure term */
     
      for(j=1;j<=ends;j++) {
	    node = E->ien[e].node[j];
	    SZZ[node] += E->TWW[E->mesh.levmax][e].node[j] * Szz;  
        }

      }

   for(node=1;node<=E->mesh.nno;node++)  {
     SZZ[node] = SZZ[node]*E->Mass[node];
     }


/*
   get_surf_stress(E,SXX,SYY,SZZ,SXY,SXZ,SZY);

    free((void *)SXX);
    free((void *)SYY);
    free((void *)SXY);
    free((void *)SXZ);
    free((void *)SZY);
*/

    return; 
}      
