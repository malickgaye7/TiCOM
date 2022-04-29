 /*   	This is where the scaling functions and grid related things are kept.

		Louis Moresi aka LUIGI   6.xii.1989                */

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


void twiddle_thumbs(yawn,scratch_groin)
     struct All_variables *yawn;
     int scratch_groin;

{ /* Do nothing, just sit back and relax.
     Take it easy for a while, maybe size
     doesn't matter after all. There, there
     that's better. Now ... */

  return; }


/*	==================================================================================
	Function to give the global shape function from the local: Assumes ORTHOGONAL MESH 
	==================================================================================      */
	
void get_global_shape_fn(E,el,GN,GNx,dOmega,x,pressure,lv)
     struct All_variables *E;
     int el;
     struct Shape_function *GN;
     struct Shape_function_dx *GNx;
     struct Shape_function_dA *dOmega;
     int pressure,lv;
     double x[3][5];
{	
  int i,j,k,d,e;
  double scale1,scale2,scale3;
  double area;
  double jacobian;
  double determinant();
  double cofactor();
 
  double dxda[4][4],cof[4][4];
 

  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];
  const int ppts=ppoints[dims];
  const int spts=spoints[dims];


       if(pressure < 2) {
           for(k=1;k<=vpts;k++) {       /* all of the vpoints */

	      for(d=1;d<=dims;d++)
		for(e=1;e<=dims;e++)
		    dxda[d][e]=0.0;

	      for(d=1;d<=dims;d++)
	         x[d][k] = 0.0; 

	      for(d=1;d<=dims;d++)
	        for(i=1;i<=ends;i++)
                  x[d][k] += E->XX[lv][d][E->IEN[lv][el].node[i]] * E->N.vpt[GNVINDEX(i,k)];

	      for(i=1;i<=ends;i++)
		 for(d=1;d<=dims;d++)
		    for(e=1;e<=dims;e++)
			dxda[d][e] += E->XX[lv][e][E->IEN[lv][el].node[i]] * E->Nx.vpt[GNVXINDEX(d-1,i,k)];   /* This is Shijie's change (d<->e) */

	      jacobian = determinant(dxda,E->mesh.nsd);  
	      dOmega->vpt[k] = jacobian*x[2][k]*x[2][k]*sin(x[1][k]);

	      for(d=1;d<=dims;d++)
		 for(e=1;e<=dims;e++)
		    cof[d][e]=cofactor(dxda,d,e,dims); 
 
	      for(j=1;j<=ends;j++)
		for(d=1;d<=dims;d++)         {
		    GNx->vpt[GNVXINDEX(d-1,j,k)] = 0.0;
		    for(e=1;e<=dims;e++)
			GNx->vpt[GNVXINDEX(d-1,j,k)] += E->Nx.vpt[GNVXINDEX(e-1,j,k)] *cof[e][d];   /* switch e and d for cof  -- Shijie's  */
		    GNx->vpt[GNVXINDEX(d-1,j,k)] /= jacobian;
	      	    }    

              x[1][k] = cos(x[1][k])/sin(x[1][k]);
              x[2][k] = 1.0/x[2][k];

	      for(j=1;j<=ends;j++)
                GNx->vpt[GNVXINDEX(0,j,k)] = GNx->vpt[GNVXINDEX(0,j,k)]*x[2][k];
	  
	      }
	    }
	
	if(pressure > 0 && pressure < 3) {
               for(k=1;k<=ppts;k++)         {   /* all of the ppoints */
		  for(d=1;d<=dims;d++)
		    for(e=1;e<=dims;e++)
			dxda[d][e]=0.0;
	          for(d=1;d<=dims;d++)
	            x[d][k] = 0.0; 

	         for(d=1;d<=dims;d++)
	           for(i=1;i<=ends;i++)
                     x[d][k] += E->XX[lv][d][E->IEN[lv][el].node[i]] * E->N.ppt[GNPINDEX(i,k)];
	      
		  for(i=1;i<=ends;i++)
		    for(d=1;d<=dims;d++)
			for(e=1;e<=dims;e++)
			    dxda[d][e] += E->XX[lv][e][E->IEN[lv][el].node[i]] * E->Nx.ppt[GNPXINDEX(d-1,i,k)];

	          jacobian = determinant(dxda,E->mesh.nsd);     
	          dOmega->ppt[k] = jacobian*x[2][k]*x[2][k]*sin(x[1][k]);

	          for(d=1;d<=dims;d++)
		     for(e=1;e<=dims;e++)
		        cof[d][e]=cofactor(dxda,d,e,E->mesh.nsd); 
	      
	          for(j=1;j<=ends;j++)
		     for(d=1;d<=dims;d++)  {
		        GNx->ppt[GNPXINDEX(d-1,j,k)]=0.0;
		        for(e=1;e<=dims;e++)
		           GNx->ppt[GNPXINDEX(d-1,j,k)] += E->Nx.ppt[GNPXINDEX(e-1,j,k)]*cof[e][d]; 
		           GNx->ppt[GNPXINDEX(d-1,j,k)] /= jacobian;
		        }

                  x[1][k] = cos(x[1][k])/sin(x[1][k]);
                  x[2][k] = 1.0/x[2][k];
	      
	          for(j=1;j<=ends;j++)
		       GNx->ppt[GNPXINDEX(0,j,k)] = GNx->ppt[GNPXINDEX(0,j,k)]*x[2][k];

	          }
	     } 


  return;
}


/*   ======================================================================
     Function to produce the appropriate one-dimensional global shape 
     function for a particular element edge. Referenced by the element/local
     node number.
     If nodal=0, then Gaussian quadrature is used; if nodal=1, then Lobatto 
     quadrature is used (i.e., calculated on nodal coordinates)
     ======================================================================  */
	
void get_global_1d_shape_fn(E,el,GM,dGammax,nodal)
     struct All_variables *E;
     int el,nodal;
     struct Shape_function1 *GM;
     struct Shape_function1_dA *dGammax;
{ 
  int i,k,d,e;
  int dirn,locn,node[5];
  int collapsed_dirn[2];
  double scale[4];

  double jacobian;
  double determinant();
  double cofactor();
  
  void get_neighbour_nodes();

  double dxda[4][4],cof[4][4];
  
  for(locn=0;locn<=1;locn++) /* top/bottom, front/back, left/right */ 
	for(dirn=1;dirn<=E->mesh.nsd;dirn++)     {
	    get_neighbour_nodes(node,dirn,locn);

	    collapsed_dirn[0]=1; collapsed_dirn[1]=3;  /* includes 3d  */
	    if (collapsed_dirn[0]==(dirn)) collapsed_dirn[0]++;
	    else if (collapsed_dirn[1]==(dirn)) collapsed_dirn[1]--;

	    for(k=1;k<=onedvpoints[E->mesh.nsd];k++)  { /* all of the vpoints*/
	        for(d=1;d<=E->mesh.nsd-1;d++)
		  for(e=1;e<=E->mesh.nsd-1;e++)
		    dxda[d][e]=0.0;
	  
		if(nodal==0)  /* Gaussian */
		   for(i=1;i<=onedvpoints[E->mesh.nsd];i++)      /* nodes */
		     for(d=1;d<=E->mesh.nsd-1;d++)
		       for(e=1;e<=E->mesh.nsd-1;e++)    {
		          dxda[d][e] += E->X[collapsed_dirn[e-1]][E->ien[el].node[node[i]]]*E->Mx.vpt[GMVXINDEX(d-1,i,k)];      
		          }
		else if(nodal==1)  /* Lobatto */
		   for(i=1;i<=onedvpoints[E->mesh.nsd];i++)      /* nodes */
		     for(d=1;d<=E->mesh.nsd-1;d++)
		       for(e=1;e<=E->mesh.nsd-1;e++)    {
		          dxda[d][e] += E->X[collapsed_dirn[e-1]][E->ien[el].node[node[i]]]*E->Lx.vpt[GMVXINDEX(d-1,i,k)];      
		          }

		jacobian = determinant(dxda,E->mesh.nsd-1); 
		dGammax->vpt[GMVGAMMA(dirn-1+E->mesh.nsd*locn,k)] = jacobian;
	        }
	  }
  
  return;
}
  
void get_neighbour_nodes(node,dirn,locn)
     int node[5];
     int dirn,locn;  /* dirn is normal to the surface, loc is front/back/top/bottom/left/right */

{ int a;  /* reference node, then use same (proven) scheme as 1d integration */

  switch(dirn)
    { case 1:    /* x vector normal */
	a = loc[1].node_nebrs[0][locn]; 
	node[1] = loc[loc[a].node_nebrs[2][0]].node_nebrs[1][0];
	node[2] = loc[loc[a].node_nebrs[2][0]].node_nebrs[1][1];
	node[4] = loc[loc[a].node_nebrs[2][1]].node_nebrs[1][0];
	node[3] = loc[loc[a].node_nebrs[2][1]].node_nebrs[1][1];
	break;
  
      case 2:    /* z vector normal */
	a = loc[1].node_nebrs[1][locn]; 
	node[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][0];
	node[2] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][0];
	node[4] = loc[loc[a].node_nebrs[0][0]].node_nebrs[2][1];
	node[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[2][1];
	break;
	
      case 3:    /* y vector normal */
	a = loc[1].node_nebrs[2][locn]; 
	node[1] = loc[loc[a].node_nebrs[0][0]].node_nebrs[1][0];
	node[2] = loc[loc[a].node_nebrs[0][0]].node_nebrs[1][1];
	node[4] = loc[loc[a].node_nebrs[0][1]].node_nebrs[1][0];
	node[3] = loc[loc[a].node_nebrs[0][1]].node_nebrs[1][1];
	break;

      }

return;
}


/*  ==========================================
    construct the lumped mass matrix. The full
    matrix is the FE integration of the density 
    field. The lumped version is the diagonal
    matrix obtained by letting the shape function
    Na be delta(a,b)
    ========================================== */

void mass_matrix(E)
     struct All_variables *E;

{ int lv,node,el,i,nint,e,le,n[9];
  void get_global_shape_fn();
  double area,centre[4],temp[9],tempa,dx1,dx2,dx3,xk[3][5];
  float xlowmean,normlow,normhigh,xhighmean;
  struct Shape_function GN;
  struct Shape_function_dA dOmega;
  struct Shape_function_dx GNx;

  const int ppts=ppoints[E->mesh.nsd];
  const int vpts=vpoints[E->mesh.nsd];

  /* ECO .size can also be defined here */


  for (lv=E->mesh.levmin;lv<=E->mesh.levmax;lv++)   {
    for(node=1;node<=E->mesh.NNO[lv];node++)
      E->MASS[lv][node] = 0.0;

    for(e=1;e<=E->mesh.NEL[lv];e++)  {

      for(node=1;node<=enodes[E->mesh.nsd];node++)
        n[node] = E->IEN[lv][e].node[node];

      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,xk,0,lv);

      area = centre[1] = centre[2] = centre[3] = 0.0;

      for(i=1;i<=E->mesh.nsd;i++)  {
        for(node=1;node<=enodes[E->mesh.nsd];node++)
	      centre[i] += E->XX[lv][i][E->IEN[lv][e].node[node]];

    	E->ECO[lv][e].centre[i] = centre[i]/enodes[E->mesh.nsd];
        }     /* end loop for dof */

      if (3==E->mesh.nsd)   {
        dx1 = 0.25*(E->XX[lv][1][n[3]]+E->XX[lv][1][n[4]]
                   +E->XX[lv][1][n[7]]+E->XX[lv][1][n[8]]
                   -E->XX[lv][1][n[1]]-E->XX[lv][1][n[2]]
                   -E->XX[lv][1][n[5]]-E->XX[lv][1][n[6]]);
        dx2 = 0.25*(E->XX[lv][2][n[3]]+E->XX[lv][2][n[4]]
                   +E->XX[lv][2][n[7]]+E->XX[lv][2][n[8]]
                   -E->XX[lv][2][n[1]]-E->XX[lv][2][n[2]]
                   -E->XX[lv][2][n[5]]-E->XX[lv][3][n[6]]);
        dx3 = 0.25*(E->XX[lv][3][n[3]]+E->XX[lv][3][n[4]]
                   +E->XX[lv][3][n[7]]+E->XX[lv][3][n[8]]
                   -E->XX[lv][3][n[1]]-E->XX[lv][3][n[2]]
                   -E->XX[lv][3][n[5]]-E->XX[lv][3][n[6]]);
        E->ECO[lv][e].size[1] = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3);

        dx1 = 0.25*(E->XX[lv][1][n[2]]+E->XX[lv][1][n[3]]
                   +E->XX[lv][1][n[6]]+E->XX[lv][1][n[7]]
                   -E->XX[lv][1][n[1]]-E->XX[lv][1][n[4]]
                   -E->XX[lv][1][n[5]]-E->XX[lv][1][n[8]]);
        dx2 = 0.25*(E->XX[lv][2][n[2]]+E->XX[lv][2][n[3]]
                   +E->XX[lv][2][n[6]]+E->XX[lv][2][n[7]]
                   -E->XX[lv][2][n[1]]-E->XX[lv][2][n[4]]
                   -E->XX[lv][2][n[5]]-E->XX[lv][3][n[8]]);
        dx3 = 0.25*(E->XX[lv][3][n[2]]+E->XX[lv][3][n[3]]
                   +E->XX[lv][3][n[6]]+E->XX[lv][3][n[7]]
                   -E->XX[lv][3][n[1]]-E->XX[lv][3][n[4]]
                   -E->XX[lv][3][n[5]]-E->XX[lv][3][n[8]]);
        E->ECO[lv][e].size[2] = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3);

        dx1 = 0.25*(E->XX[lv][1][n[5]]+E->XX[lv][1][n[6]]
                   +E->XX[lv][1][n[7]]+E->XX[lv][1][n[8]]
                   -E->XX[lv][1][n[1]]-E->XX[lv][1][n[2]]
                   -E->XX[lv][1][n[3]]-E->XX[lv][1][n[4]]);
        dx2 = 0.25*(E->XX[lv][2][n[5]]+E->XX[lv][2][n[6]]
                   +E->XX[lv][2][n[7]]+E->XX[lv][2][n[8]]
                   -E->XX[lv][2][n[1]]-E->XX[lv][2][n[2]]
                   -E->XX[lv][2][n[3]]-E->XX[lv][3][n[4]]);
        dx3 = 0.25*(E->XX[lv][3][n[5]]+E->XX[lv][3][n[6]]
                   +E->XX[lv][3][n[7]]+E->XX[lv][3][n[8]]
                   -E->XX[lv][3][n[1]]-E->XX[lv][3][n[2]]
                   -E->XX[lv][3][n[3]]-E->XX[lv][3][n[4]]);
        E->ECO[lv][e].size[3] = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3);

        }
      else if (2==E->mesh.nsd)   {
        dx1 = 0.5*(E->XX[lv][1][n[3]]+E->XX[lv][1][n[4]]
                  -E->XX[lv][1][n[1]]-E->XX[lv][1][n[2]]);
        dx2 = 0.5*(E->XX[lv][2][n[3]]+E->XX[lv][2][n[4]]
                  -E->XX[lv][2][n[1]]-E->XX[lv][2][n[2]]);
        E->ECO[lv][e].size[1] = sqrt(dx1*dx1 + dx2*dx2);

        dx1 = 0.5*(E->XX[lv][1][n[2]]+E->XX[lv][1][n[3]]
                  -E->XX[lv][1][n[1]]-E->XX[lv][1][n[4]]);
        dx2 = 0.5*(E->XX[lv][2][n[2]]+E->XX[lv][2][n[3]]
                  -E->XX[lv][2][n[1]]-E->XX[lv][2][n[4]]);
        E->ECO[lv][e].size[2] = sqrt(dx1*dx1 + dx2*dx2);
        }

      for(nint=1;nint<=vpts;nint++)
        area += g_point[nint].weight[E->mesh.nsd-1] * dOmega.vpt[nint];
      E->ECO[lv][e].area = area;

      for(node=1;node<=enodes[E->mesh.nsd];node++)  {
        temp[node] = 0.0;
        for(nint=1;nint<=vpts;nint++)
          temp[node] += dOmega.vpt[nint]*g_point[nint].weight[E->mesh.nsd-1]
                       *E->N.vpt[GNVINDEX(node,nint)];       /* int Na dV */
        }

      for(node=1;node<=enodes[E->mesh.nsd];node++)
         E->MASS[lv][E->IEN[lv][e].node[node]] += temp[node];

        for(node=1;node<=enodes[E->mesh.nsd];node++)  {
           E->TWW[lv][e].node[node] = temp[node];
           }

      } /* over element */
    for(node=1;node<=E->mesh.NNO[lv];node++)  {
      E->MASS[lv][node] = 1.0/E->MASS[lv][node];
      }


    }

	E->Total.vol = 0.0;
  for(e=1;e<E->mesh.nel;e++)  
    E->Total.vol += E->ECO[E->mesh.levmax][e].area;
     
  fprintf(E->fp,"Total mantle volume: %.5e\n",E->Total.vol);

  return;
 }
