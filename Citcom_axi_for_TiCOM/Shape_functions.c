/*  Functions which construct the shape function values at all of the gauss
    points in the element (including the reduced quadrature points). The element in question is
    biquadratic in the velocities and therefore bilinear in the pressures. 
    
    To change elements it is necessary to change this file: Shape_functions.c,
    and the element-data header file : element_definitions.h  but it should not be
    necessary to change the main calculation/setup/solving machinery.		 */

#include <math.h>
#include "element_definitions.h"				
#include "global_defs.h"
 
/*  =======================================================
    Function creating shape_fn data in form of a structure
    =======================================================*/

void construct_shape_functions(E)
     struct All_variables *E;
{	
  double lpoly(),lpolydash();
  int i,j,k,d,dd;
  int remapj,remapk;

  /* first zero ALL entries, even those not used in 2d. */

  for(i=0;i<GNVI;i++)
    { E->N.vpt[i] = 0.0; 
      E->Nx.vpt[i] = 0.0;
      E->Nx.vpt[GNVI+i] = 0.0;
      E->Nx.vpt[2*GNVI+i] = 0.0; 
    }
  
   for(i=0;i<GNPI;i++)
    { E->N.ppt[i] = 0.0; 
      E->Nx.ppt[i] = 0.0;
      E->Nx.ppt[GNPI+i] = 0.0;
      E->Nx.ppt[2*GNPI+i] = 0.0; 
    }
  
  for(i=0;i<GN1VI;i++)
    { E->M.vpt[i] = 0.0; 
      E->Mx.vpt[i] = 0.0;
      E->Mx.vpt[GN1VI+i] = 0.0;
    }
  
   for(i=0;i<GN1PI;i++)
    { E->M.ppt[i] = 0.0; 
      E->Mx.ppt[i] = 0.0;
      E->Mx.ppt[GN1PI+i] = 0.0;
    }
  
  for(i=0;i<GN1VI;i++)
    { E->L.vpt[i] = 0.0; 
      E->Lx.vpt[i] = 0.0;
      E->Lx.vpt[GN1VI+i] = 0.0;
    }

  for(i=1;i<=enodes[E->mesh.nsd];i++)   {
   /*  for each node  */

      for(j=1;j<=vpoints[E->mesh.nsd];j++)  { 
 
	  /* for each integration point  */
	  E->N.vpt[GNVINDEX(i,j)] = 1.0;
	  for(d=1;d<=E->mesh.nsd;d++)  {
	      E->N.vpt[GNVINDEX(i,j)] *=  
		  lpoly(bb[d-1][i],g_point[j].x[d-1]);
	  }
	  for(dd=1;dd<=E->mesh.nsd;dd++) {
	      E->Nx.vpt[GNVXINDEX(dd-1,i,j)] = lpolydash(bb[dd-1][i],g_point[j].x[dd-1]);
	      for(d=1;d<=E->mesh.nsd;d++)
		  if (d != dd)
		      E->Nx.vpt[GNVXINDEX(dd-1,i,j)] *= lpoly(bb[d-1][i],g_point[j].x[d-1]);
	  }
      } 
 
     
      for(j=1;j<=ppoints[E->mesh.nsd];j++)  {
	  /* for each p-integration point  */
	  E->N.ppt[GNPINDEX(i,j)] = 1.0;
	  for(d=1;d<=E->mesh.nsd;d++){
	      E->N.ppt[GNPINDEX(i,j)] *=  
		  lpoly(bb[d-1][i],p_point[j].x[d-1]);
	  }
	  for(dd=1;dd<=E->mesh.nsd;dd++) {
	      E->Nx.ppt[GNPXINDEX(dd-1,i,j)] = lpolydash(bb[dd-1][i],p_point[j].x[dd-1]);
	      for(d=1;d<=E->mesh.nsd;d++)
		  if (d != dd) {
		      E->Nx.ppt[GNPXINDEX(dd-1,i,j)] *= lpoly(bb[d-1][i],p_point[j].x[d-1]); 
		  }
	  }
      }
      
  }	 


  /*	1d cases ... set up M (Mx is not used in FE formulation but is used for dxdxsi &c)  */

  /* note: naughty boy wedgie, this works best for spatial ordering (z,x,y)
     so we actually need to read the z part if it's one-d and the x,z part if
     it's two. ' Comes about from the node numbering really and so we need the
     remapping of j in the 1d case. */

  for(j=1;j<=onedvpoints[E->mesh.nsd];j++)
    { remapj = ccc[E->mesh.nsd-2][j];
      for(k=1;k<=onedvpoints[E->mesh.nsd];k++)     
	{ remapk = ccc[E->mesh.nsd-2][k];
	  E->M.vpt[GMVINDEX(j,k)] = 1.0;
	  E->L.vpt[GMVINDEX(j,k)] = 1.0;
	  for(d=1;d<=E->mesh.nsd-1;d++) {
	    E->M.vpt[GMVINDEX(j,k)] *= lpoly(cc[d-1][remapj],g_1d[k].x[d-1]);
	    E->L.vpt[GMVINDEX(j,k)] *= lpoly(cc[d-1][remapj],l_1d[k].x[d-1]);
	    }
	  for(dd=1;dd<=E->mesh.nsd-1;dd++) {
	      E->Mx.vpt[GMVXINDEX(dd-1,j,k)] = lpolydash(cc[dd-1][remapj],g_1d[k].x[dd-1]); 
	      E->Lx.vpt[GMVXINDEX(dd-1,j,k)] = lpolydash(cc[dd-1][remapj],l_1d[k].x[dd-1]); 
	      for(d=1;d<=E->mesh.nsd-1;d++)
		if (d != dd) {
		  E->Mx.vpt[GMVXINDEX(dd-1,j,k)] *= lpoly(cc[d-1][remapj],g_1d[k].x[d-1]); 
		  E->Lx.vpt[GMVXINDEX(dd-1,j,k)] *= lpoly(cc[d-1][remapj],l_1d[k].x[d-1]); 
		  }
	      }
	  }
      
    }
    	
  return; }
		
double lpoly(p,y)
     int p;	   /*   selects lagrange polynomial , 1d: node p */
     double y;  /*   coordinate in given direction to evaluate poly */
{	
  double value;
  
  switch (p)
    {
    case 1:
      value =0.5 * (1-y) ;
      break;
    case 2:
      value =0.5 * (1+y) ;
      break;
    default:
      value = 0.0;
    }

  return(value);
}
	
double lpolydash(p,y)
     int p;
     double y;
{	
  double value;
  switch (p)
    {
    case 1:
      value = -0.5 ;
      break;
    case 2:
      value =  0.5 ;
      break;
    default:
      value = 0.0;
    }

  return(value);	}










