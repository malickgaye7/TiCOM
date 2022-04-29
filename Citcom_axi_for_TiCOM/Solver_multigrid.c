#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>

static double weightIJ[5][5] =
        {       { 0.0,    0.0,    0.0,    0.0,    0.0 },
		{ 0.0, 0.5625, 0.1875, 0.0625, 0.1875 },
		{ 0.0, 0.1875, 0.5625, 0.1875, 0.0625 },
		{ 0.0, 0.0625, 0.1875, 0.5625, 0.1875 },
		{ 0.0, 0.1875, 0.0625, 0.1875, 0.5625 } };


void set_mg_defaults(E)
     struct All_variables *E;
{ void assemble_forces_iterative();
  void solve_constrained_flow_iterative();
  void mg_allocate_vars();

  E->solver_allocate_vars = mg_allocate_vars;
  E->build_forcing_term = assemble_forces_iterative;
  E->solve_stokes_problem = solve_constrained_flow_iterative;

  E->control.mg_cycle = 1;
 
return;
}

void mg_allocate_vars(E)
     struct All_variables *E;
{  
  return;

}

void project_vector_a(E,start_lev,AU,AD,ic)
     struct All_variables *E;
     int start_lev,ic;
     double *AU,*AD;  /* data on upper/lower mesh  */
{
    int i,j;
    int el,node,node1,e1;
    int eqn1,eqn_minus1;
    int eqn2,eqn_minus2;
    int eqn3,eqn_minus3;
    double amplifier,average1,average2,average3,w;
    float CPU_time(),time;
 
    const int sl_minus = start_lev-1;
    const int neq_minus=E->mesh.NEQ[start_lev-1];
    const int nno_minus=E->mesh.NNO[start_lev-1];
    const int nels_minus=E->mesh.NEL[start_lev-1];
    const int  dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];
    const double weight=(double) 1.0/ends;
  

    /* on the lower level the average value of data in upper level
       ELEMENTS are slid across to the nodes.  */
 
   for(i=0;i<=neq_minus;i++)
      AD[i] = 0.0;

    if(3==dims)
	for(el=1;el<=nels_minus;el++) {
	    e1 = E->EL[sl_minus][el].sub[1];
            node1=E->IEN[start_lev][e1].node[7];
	    for(i=1;i<=ENODES3D;i++)                  {
		node= E->IEN[sl_minus][el].node[i];

		AD[E->ID[sl_minus][node].doff[1]] += AU[E->ID[start_lev][node1].doff[1]]; 
		AD[E->ID[sl_minus][node].doff[2]] += AU[E->ID[start_lev][node1].doff[2]]; 
		AD[E->ID[sl_minus][node].doff[3]] += AU[E->ID[start_lev][node1].doff[3]]; 
                }
	 
	 }
    else  
	for(el=1;el<=nels_minus;el++)   {

	  for(i=1;i<=ENODES2D;i++)                 {
	    e1 = E->EL[sl_minus][el].sub[i];
            average1 = average2 = 0.0;
	    for(j=1;j<=ENODES2D;j++)                 {
              node1 = E->IEN[start_lev][e1].node[j];
              average1 += AU[E->ID[start_lev][node1].doff[1]];
              average2 += AU[E->ID[start_lev][node1].doff[2]];
              }
            average1 *= weight;
            average2 *= weight;

	    node= E->IEN[sl_minus][el].node[i];

	    AD[E->ID[sl_minus][node].doff[1]] += average1; 
	    AD[E->ID[sl_minus][node].doff[2]] += average2; 
            }
	 
	 }
return;
}

void project_vector_diag(E,start_lev,AU,AD,ic)
     struct All_variables *E;
     int start_lev,ic;
     double *AU,*AD;  /* data on upper/lower mesh  */
{
    int i,j;
    int el,node,node1,e1;
    int eqn1,eqn_minus1;
    int eqn2,eqn_minus2;
    int eqn3,eqn_minus3;
    double amplifier,average1,average2,average3,w;
    float CPU_time(),time;
 
    const int sl_minus = start_lev-1;
    const int neq_minus=E->mesh.NEQ[start_lev-1];
    const int nno_minus=E->mesh.NNO[start_lev-1];
    const int nels_minus=E->mesh.NEL[start_lev-1];
    const int  dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];
    const double weight=(double) 1.0/ends;
  

    /* on the lower level the average value of data in upper level
       ELEMENTS are slid across to the nodes.  */
 
   for(i=0;i<=neq_minus;i++)
      AD[i] = 0.0;

    if(3==dims)
	for(el=1;el<=nels_minus;el++) {
	    e1 = E->EL[sl_minus][el].sub[1];
            node1=E->IEN[start_lev][e1].node[7];
	    for(i=1;i<=ENODES3D;i++)                  {
		node= E->IEN[sl_minus][el].node[i];

		AD[E->ID[sl_minus][node].doff[1]] += AU[E->ID[start_lev][node1].doff[1]]; 
		AD[E->ID[sl_minus][node].doff[2]] += AU[E->ID[start_lev][node1].doff[2]]; 
		AD[E->ID[sl_minus][node].doff[3]] += AU[E->ID[start_lev][node1].doff[3]]; 
                }
	 
	 }
    else  
	for(el=1;el<=nels_minus;el++) {
	    e1 = E->EL[sl_minus][el].sub[1];
            node1=E->IEN[start_lev][e1].node[3];
	    for(i=1;i<=ENODES2D;i++)                 {
		node= E->IEN[sl_minus][el].node[i];

		AD[E->ID[sl_minus][node].doff[1]] += AU[E->ID[start_lev][node1].doff[1]]; 
		AD[E->ID[sl_minus][node].doff[2]] += AU[E->ID[start_lev][node1].doff[2]]; 
                }
	 
	 }
return;
}


void project_vector(E,start_lev,AU,AD,ic)
     struct All_variables *E;
     int start_lev,ic;
     double *AU,*AD;  /* data on upper/lower mesh  */
{
    int i,j;
    int el,node,e1;
    int eqn1,eqn_minus1;
    int eqn2,eqn_minus2;
    int eqn3,eqn_minus3;
    double amplifier,average1,average2,average3,w;
    float CPU_time(),time;
 
    const int sl_minus = start_lev-1;
    const int neq_minus=E->mesh.NEQ[start_lev-1];
    const int nno_minus=E->mesh.NNO[start_lev-1];
    const int nels_minus=E->mesh.NEL[start_lev-1];
    const int  dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];
    const double weight=(double) 1.0/ends;
  

    /* on the lower level the average value of data in upper level
       ELEMENTS are slid across to the nodes.  */
 
   for(i=0;i<=neq_minus;i++)
      AD[i] = 0.0;

    if(3==dims)
	for(el=1;el<=nels_minus;el++)
	    for(i=1;i<=ENODES3D;i++) {
		average1=average2=average3=0.0;
		e1 = E->EL[sl_minus][el].sub[i];
		for(j=1;j<=ENODES3D;j++) {
		    node=E->IEN[start_lev][e1].node[j];
		    average1 += AU[E->ID[start_lev][node].doff[1]];
		    average2 += AU[E->ID[start_lev][node].doff[2]];
		    average3 += AU[E->ID[start_lev][node].doff[3]];
	      
		}     
		node= E->IEN[sl_minus][el].node[i];

		AD[E->ID[sl_minus][node].doff[1]] += E->TWW[sl_minus][el].node[i] * average1; 
		AD[E->ID[sl_minus][node].doff[2]] += E->TWW[sl_minus][el].node[i] * average2; 
	 	AD[E->ID[sl_minus][node].doff[3]] += E->TWW[sl_minus][el].node[i] * average3; 
	 
	 }
    else  
	for(el=1;el<=nels_minus;el++)
	    for(i=1;i<=ENODES2D;i++)	{ 
		average1=average2=0.0;
		e1 = E->EL[sl_minus][el].sub[i];
		for(j=1;j<=ENODES2D;j++)    {
		    node =E->IEN[start_lev][e1].node[j];
		    average1 += AU[E->ID[start_lev][node].doff[1]];
		    average2 += AU[E->ID[start_lev][node].doff[2]];
		}    	
		node=E->IEN[sl_minus][el].node[i];
		AD[E->ID[sl_minus][node].doff[1]] += E->TWW[sl_minus][el].node[i] * average1; 
		AD[E->ID[sl_minus][node].doff[2]] += E->TWW[sl_minus][el].node[i] * average2; 
	    }

   if(3==dims)
     for(i=1;i<=nno_minus;i++)  {
       AD[E->ID[sl_minus][i].doff[1]] = AD[E->ID[sl_minus][i].doff[1]] * E->MASS[sl_minus][i];
       AD[E->ID[sl_minus][i].doff[2]] = AD[E->ID[sl_minus][i].doff[2]] * E->MASS[sl_minus][i];
       AD[E->ID[sl_minus][i].doff[3]] = AD[E->ID[sl_minus][i].doff[3]] * E->MASS[sl_minus][i];
       }
   else
     for(i=1;i<=nno_minus;i++)  {
       AD[E->ID[sl_minus][i].doff[1]] = AD[E->ID[sl_minus][i].doff[1]] * E->MASS[sl_minus][i];
       AD[E->ID[sl_minus][i].doff[2]] = AD[E->ID[sl_minus][i].doff[2]] * E->MASS[sl_minus][i];
       }
    
/* Thats's all */

return;  }

void project_vector1(E,start_lev,AU,AD,ic)

     struct All_variables *E;
     int start_lev,ic;
     double *AU,*AD;  /* data on upper/lower mesh  */
{
    int i,j;
    int el,node,e1;
    int eqn1,eqn_minus1;
    int eqn2,eqn_minus2;
    int eqn3,eqn_minus3;
    double amplifier,average1,average2,average3,w;
    float CPU_time(),time;
 
    const int sl_minus = start_lev-1;
    const int neq_minus=E->mesh.NEQ[start_lev-1];
    const int nno_minus=E->mesh.NNO[start_lev-1];
    const int nels_minus=E->mesh.NEL[start_lev-1];
    const int  dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];
    const double weight=(double) 1.0/ends;
  

    /* on the lower level the average value of data in upper level
       ELEMENTS are slid across to the nodes.  */
 
   for(i=0;i<=neq_minus;i++)
      AD[i] = 0.0;

    if(3==dims)
	for(el=1;el<=nels_minus;el++)
	    for(i=1;i<=ENODES3D;i++) {
		average1=average2=average3=0.0;
		e1 = E->EL[sl_minus][el].sub[i];
		for(j=1;j<=ENODES3D;j++) {
		    node=E->IEN[start_lev][e1].node[j];
		    average1 += AU[E->ID[start_lev][node].doff[1]];
		    average2 += AU[E->ID[start_lev][node].doff[2]];
		    average3 += AU[E->ID[start_lev][node].doff[3]];
	      
		}     
		node= E->IEN[sl_minus][el].node[i];
		w=weight * E->TWW[sl_minus][el].node[i]; 

		AD[E->ID[sl_minus][node].doff[1]] += w * average1; 
		AD[E->ID[sl_minus][node].doff[2]] += w * average2; 
	 	AD[E->ID[sl_minus][node].doff[3]] += w * average3; 
	 
	 }
    else  
	for(el=1;el<=nels_minus;el++)
	    for(i=1;i<=ENODES2D;i++)	{ 
		average1=average2=0.0;
		e1 = E->EL[sl_minus][el].sub[i];
		for(j=1;j<=ENODES2D;j++)    {
		    node =E->IEN[start_lev][e1].node[j];
		    average1 += AU[E->ID[start_lev][node].doff[1]];
		    average2 += AU[E->ID[start_lev][node].doff[2]];
		}    	
		node=E->IEN[sl_minus][el].node[i];
		w=weight * E->TWW[sl_minus][el].node[i]; 
		AD[E->ID[sl_minus][node].doff[1]] += w * average1; 
		AD[E->ID[sl_minus][node].doff[2]] += w * average2; 
	    }

   if (ic==1) amplifier = 4.0;
   else amplifier = 1.0;

   if(3==dims)
     for(i=1;i<=nno_minus;i++)  {
       AD[E->ID[sl_minus][i].doff[1]] = AD[E->ID[sl_minus][i].doff[1]] * E->MASS[sl_minus][i] * amplifier;
       AD[E->ID[sl_minus][i].doff[2]] = AD[E->ID[sl_minus][i].doff[2]] * E->MASS[sl_minus][i] * amplifier;
       AD[E->ID[sl_minus][i].doff[3]] = AD[E->ID[sl_minus][i].doff[3]] * E->MASS[sl_minus][i] * amplifier;
       }
   else
     for(i=1;i<=nno_minus;i++)  {
       AD[E->ID[sl_minus][i].doff[1]] = AD[E->ID[sl_minus][i].doff[1]] * E->MASS[sl_minus][i] * amplifier;
       AD[E->ID[sl_minus][i].doff[2]] = AD[E->ID[sl_minus][i].doff[2]] * E->MASS[sl_minus][i] * amplifier;
       }
    
/* Thats's all */

return;  }

/* =====================================================
   Function to inject data from high to low grid (i.e.
   just dropping values not at shared grid points.
   ===================================================== */

void inject(E,start_lev,AU,AD)

     struct All_variables *E;
     int start_lev;
     double *AU,*AD;  /* data on upper/lower mesh  */

{   
    int i;
    int el,node_coarse,node_fine;
    int sl_minus;
    int eqn,eqn_coarse;
  
    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];

    if(start_lev == E->mesh.levmin)   {
	fprintf(E->fp,"Warning, attempting to project below lowest level\n");
	return;
    }
    sl_minus = start_lev-1;

    for(el=1;el<=E->mesh.NEL[sl_minus];el++)
	for(i=1;i<=ends;i++) {
	    node_coarse = E->IEN[sl_minus][el].node[i];
	    node_fine=E->IEN[start_lev][E->EL[sl_minus][el].sub[i]].node[i];
	   
	    eqn_coarse = E->ID[sl_minus][node_coarse].doff[1];
	    eqn = E->ID[start_lev][node_fine].doff[1];
	    AD[eqn_coarse] = AU[eqn]; 
	   
	    eqn_coarse = E->ID[sl_minus][node_coarse].doff[2];
	    eqn = E->ID[start_lev][node_fine].doff[2];
	    AD[eqn_coarse] = AU[eqn]; 
	   
	    if(3==dims)  {
		eqn_coarse = E->ID[sl_minus][node_coarse].doff[3];
		eqn = E->ID[start_lev][node_fine].doff[3];
		AD[eqn_coarse] = AU[eqn]; 
	    }
	} 
    return; 
}

/* =====================================================
   Function to inject data from high to low grid (i.e.
   just dropping values not at shared grid points.
   ===================================================== */

void un_inject_vector(E,start_lev,AD,AU)

     struct All_variables *E;
     int start_lev;
     double *AU,*AD;  /* data on upper/lower mesh  */
{  
    int i;
    int el,node,node_plus;
    int eqn1,eqn_plus1;
    int eqn2,eqn_plus2;
    int eqn3,eqn_plus3;
    
   
    
    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int sl_plus = start_lev+1;
    const int neq = E->mesh.NEQ[sl_plus];
    const int nels = E->mesh.NEL[start_lev];

    assert(start_lev != E->mesh.levmax  /* un_injection */);

    for(i=1;i<=neq;i++)
	AU[i]=0.0;
    
    if(3==dims)
	for(el=1;el<=nels;el++)
	    for(i=1;i<=ENODES3D;i++)  {
		node = E->IEN[start_lev][el].node[i];
		node_plus=E->IEN[sl_plus][E->EL[start_lev][el].sub[i]].node[i];
	    
		eqn1 = E->ID[start_lev][node].doff[1];
		eqn2 = E->ID[start_lev][node].doff[2];
		eqn3 = E->ID[start_lev][node].doff[3];
		eqn_plus1 = E->ID[sl_plus][node_plus].doff[1];
		eqn_plus2 = E->ID[sl_plus][node_plus].doff[2];
		eqn_plus3 = E->ID[sl_plus][node_plus].doff[3];
		AU[eqn_plus1] = AD[eqn1];      
		AU[eqn_plus2] = AD[eqn2]; 
		AU[eqn_plus3] = AD[eqn3]; 
	
	    } 
    else
	for(el=1;el<=nels;el++)
	    for(i=1;i<=ENODES2D;i++)  {
		node = E->IEN[start_lev][el].node[i];
		node_plus=E->IEN[sl_plus][E->EL[start_lev][el].sub[i]].node[i];
	    
		eqn1 = E->ID[start_lev][node].doff[1];
		eqn2 = E->ID[start_lev][node].doff[2];
		eqn_plus1 = E->ID[sl_plus][node_plus].doff[1];
		eqn_plus2 = E->ID[sl_plus][node_plus].doff[2];
		AU[eqn_plus1] = AD[eqn1];
		AU[eqn_plus2] = AD[eqn2]; 
	    } 	

    return; 
}

void inject_scalar(E,start_lev,AU,AD)
     struct All_variables *E;
     int start_lev;
     float *AU,*AD;  /* data on upper/lower mesh  */

{
    int i,m,el,node_coarse,node_fine,sl_minus,eqn,eqn_coarse;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];

    if(start_lev == E->mesh.levmin)   {
        fprintf(E->fp,"Warning, attempting to project below lowest level\n");
        return;
    }

    sl_minus = start_lev-1;

      for(el=1;el<=E->mesh.NEL[sl_minus];el++)   {
        for(i=1;i<=ends;i++)       {
          node_coarse = E->IEN[sl_minus][el].node[i];
          node_fine=E->IEN[start_lev][E->EL[sl_minus][el].sub[i]].node[i];
          AD[node_coarse] = AU[node_fine];
          }
        }

    return;
}


void inject_node_fvector(E,start_lev,AU,AD)

     struct All_variables *E;
     int start_lev;
     float **AU,**AD;  /* data on upper/lower mesh  */
{
    int i;
    int el,ex,ey,ez,d;
    int node,node_minus;
    
    const int sl_minus = start_lev-1;
    const int elx=E->mesh.ELX[sl_minus];
    const int elz=E->mesh.ELZ[sl_minus];
    const int ely=E->mesh.ELY[sl_minus];
    const int dims=E->mesh.nsd;
    
    assert(start_lev != E->mesh.levmin );
   
  	for(ey=1;ey<=ely;ey++)
	    for(ez=1;ez<=elz;ez++) 
		for(ex=1;ex<=elx;ex++) {
		    el=ez+(ex-1)*elz+(ey-1)*elz*elx;
		    for(i=1;i<=enodes[dims];i++) {
			node_minus = E->IEN[sl_minus][el].node[i];
			node = E->IEN[start_lev][E->EL[sl_minus][el].sub[i]].node[i];
			
			if(E->mesh.periodic_x && elx==ex && 0==loc[i].plus[0]) {
			    node_minus += (E->mesh.NOX[sl_minus]-1)*E->mesh.NOZ[sl_minus];
			    node += (E->mesh.NOX[start_lev]-1)*E->mesh.NOZ[start_lev];
			}
			AD[1][node_minus] = AU[1][node]; 
			AD[2][node_minus] = AU[2][node]; 
			if(3==dims)
			    AD[3][node_minus] = AU[3][node]; 
		    } 
		}
    return; 
}

/* =======================================================================================
   Interpolation from coarse grid to fine. See the appology attached to project() if you get
   stressed out by node based assumptions. If it makes you feel any better, I don't like
   it much either.
   ======================================================================================= */


void interp_vector(E,start_lev,AD,AU)

    struct All_variables *E;
     int start_lev;
     double *AD,*AU;  /* data on upper/lower mesh  */
{
    void un_inject_vector();
    int i,j,k;
    float x1,x2;
    float n1,n2;
    int node0,node1,node2;
    int eqn0,eqn1,eqn2;

    const int level = start_lev + 1;
    const int dims =E->mesh.nsd;
    const int ends= enodes[dims];
    
    const int nox = E->mesh.NOX[level];
    const int noz = E->mesh.NOZ[level];
    const int noy = E->mesh.NOY[level];
    const int high_eqn = E->mesh.NEQ[level];
 
    if (start_lev==E->mesh.levmax) return;

    un_inject_vector(E,start_lev,AD,AU); /* transfer information from lower level */

    for(k=1;k<=noy;k+=2)          /* Fill in gaps in x direction */
      for(j=1;j<=noz;j+=2)
	  for(i=2;i<nox;i+=2)  {
	      node0 = j + (i-1)*noz + (k-1)*nox*noz; /* this node */
	      node1 = node0 - noz;
	      node2 = node0 + noz;

  	      x1=E->Interp[level][1][node1];
	      x2=E->Interp[level][1][node2];
		
	      n1=x2/(x1+x2);
	      n2=x1/(x1+x2);
	      
	      /* now for each direction */
	      
	      eqn0=E->ID[level][node0].doff[1];
	      eqn1=E->ID[level][node1].doff[1];
	      eqn2=E->ID[level][node2].doff[1]; 
	      AU[eqn0] = n1*AU[eqn1]+n2*AU[eqn2];

	      eqn0=E->ID[level][node0].doff[2];
	      eqn1=E->ID[level][node1].doff[2];
	      eqn2=E->ID[level][node2].doff[2];
	      AU[eqn0] = n1* AU[eqn1]+n2*AU[eqn2];

	      if(3==dims)  {
		  eqn0=E->ID[level][node0].doff[3];
		  eqn1=E->ID[level][node1].doff[3];
		  eqn2=E->ID[level][node2].doff[3];	       
		  AU[eqn0] = n1*AU[eqn1]+n2*AU[eqn2];
	      }
	      
	}


    for(k=1;k<=noy;k+=2)          /* Fill in gaps in z direction */
	for(i=1;i<=nox;i++)
	    for(j=2;j<noz;j+=2)	{ 
		node0 = j + (i-1)*noz + (k-1)*nox*noz; /* this node */
		node1 = node0 - 1;
		node2 = node0 + 1;

  	        x1=E->Interp[level][2][node1];
	        x2=E->Interp[level][2][node2];
	
		n1=x2/(x1+x2);
		n2=x1/(x1+x2);
 
		eqn0=E->ID[level][node0].doff[1];
		eqn1=E->ID[level][node1].doff[1];
		eqn2=E->ID[level][node2].doff[1];  
		AU[eqn0] = n1*AU[eqn1]+n2*AU[eqn2];

		eqn0=E->ID[level][node0].doff[2];
		eqn1=E->ID[level][node1].doff[2];
		eqn2=E->ID[level][node2].doff[2];
		AU[eqn0] = n1*AU[eqn1]+n2*AU[eqn2];

		if(3==dims)  {
		    eqn0=E->ID[level][node0].doff[3];
		    eqn1=E->ID[level][node1].doff[3];
		    eqn2=E->ID[level][node2].doff[3];	       
		    AU[eqn0] = n1*AU[eqn1]+n2*AU[eqn2];
		}
	    }

          
    if(3==dims)
	for(i=1;i<=nox;i++)   /* Fill in gaps in y direction */
	    for(j=1;j<=noz;j++)
		for(k=2;k<noy;k+=2)   {
		    node0 = j + (i-1)*noz + (k-1)*nox*noz; /* this node */
		    node1 = node0 - nox*noz;
		    node2 = node0 + nox*noz;

	            x1=E->Interp[level][3][node1];
	            x2=E->Interp[level][3][node2];
	  		   
		    n1=x2/(x1+x2);
		    n2=x1/(x1+x2);
  
		    eqn0=E->ID[level][node0].doff[1];
		    eqn1=E->ID[level][node1].doff[1];
		    eqn2=E->ID[level][node2].doff[1];  
		    AU[eqn0] = n1*AU[eqn1]+n2*AU[eqn2];
	    
		    eqn0=E->ID[level][node0].doff[2];
		    eqn1=E->ID[level][node1].doff[2];
		    eqn2=E->ID[level][node2].doff[2];
		    AU[eqn0] = n1*AU[eqn1]+n2*AU[eqn2];

		    eqn0=E->ID[level][node0].doff[3];
		    eqn1=E->ID[level][node1].doff[3];
		    eqn2=E->ID[level][node2].doff[3];	       
		    AU[eqn0] = n1*AU[eqn1]+n2*AU[eqn2];	    
	}
  return;

}


/*  ==============================================
    function to project viscosity down to all the 
    levels in the problem. (no gaps for vbcs)
    ==============================================  */

void project_viscosity(E)
     struct All_variables *E;

{ 
    int lv,i,j,k,el,sl_minus;
    void inject_scalar();
    void project_scalar();
    void project_scalar_e();
    void inject_scalar_e();
    void visc_from_gint_to_nodes();
    void visc_from_nodes_to_gint();
    void visc_from_gint_to_ele();
    void visc_from_ele_to_gint();
    void build_interp_stencil();

    const int nsd=E->mesh.nsd;
    const int vpts=vpoints[nsd];
   
    float *viscU,*viscD;

    viscU=(float *)malloc((1+vpts*E->mesh.NEL[E->mesh.levmax  ])*sizeof(float));
    viscD=(float *)malloc((1+vpts*E->mesh.NEL[E->mesh.levmax-1])*sizeof(float));
/*
    visc_from_gint_to_nodes(E,E->EVI[E->mesh.levmax],viscU,E->mesh.levmax);
    visc_from_nodes_to_gint(E,viscU,E->EVI[E->mesh.levmax],E->mesh.levmax);
*/
  for(lv=E->mesh.levmax;lv>E->mesh.levmin;lv--){
    sl_minus = lv -1;

    if (E->viscosity.smooth_cycles==0)  {
      visc_from_gint_to_nodes(E,E->EVI[lv],viscU,lv);
      inject_scalar(E,lv,viscU,viscD);
      visc_from_nodes_to_gint(E,viscD,E->EVI[sl_minus],sl_minus);
      }
    else if (E->viscosity.smooth_cycles==1)  {
      visc_from_gint_to_nodes(E,E->EVI[lv],viscU,lv);
      project_scalar(E,lv,viscU,viscD);
      visc_from_nodes_to_gint(E,viscD,E->EVI[sl_minus],sl_minus);
      }
    else if (E->viscosity.smooth_cycles==2)   {
      visc_from_gint_to_ele(E,E->EVI[lv],viscU,lv);
      inject_scalar_e(E,lv,viscU,E->EVI[sl_minus]);
      }
    else if (E->viscosity.smooth_cycles==3)   {
      visc_from_gint_to_ele(E,E->EVI[lv],viscU,lv);
      project_scalar_e(E,lv,viscU,viscD);
      visc_from_ele_to_gint(E,viscD,E->EVI[sl_minus],sl_minus);
      }

    visc_from_gint_to_nodes(E,E->EVI[lv],viscU,lv);
    build_interp_stencil(E,viscU,lv);

    }

    free((void *)viscU);
    free((void *)viscD);
    

    return;  
}

void build_interp_stencil(E,viscU,lev)
 struct All_variables *E;
 float *viscU;
 int lev;
{
    const int dims =E->mesh.nsd;
    const int ends= enodes[dims];

    const int nox = E->mesh.NOX[lev];
    const int noz = E->mesh.NOZ[lev];
    const int noy = E->mesh.NOY[lev];

 float x1,x2;
 int i,j,k,node0,node1,node2;

  for(k=1;k<=noy;k+=2)          /* Fill in gaps in x direction */
    for(j=1;j<=noz;j+=2)
       for(i=2;i<nox;i+=2)  {
           node0 = j + (i-1)*noz + (k-1)*nox*noz; /* this node */
           node1 = node0 - noz;
           node2 = node0 + noz;
           x1=E->ECO[lev][E->NEI[lev].element[ends*(node1-1)]].size[1];
           x2=E->ECO[lev][E->NEI[lev].element[(node2-1)*ends]].size[1];

           E->Interp[lev][1][node1] = x1;
           E->Interp[lev][1][node2] = x2;
           }

  for(k=1;k<=noy;k+=2)          /* Fill in gaps in z direction */
    for(i=1;i<=nox;i++)
       for(j=2;j<noz;j+=2) {
           node0 = j + (i-1)*noz + (k-1)*nox*noz; /* this node */
           node1 = node0 - 1;
           node2 = node0 + 1;
           x1=E->ECO[lev][E->NEI[lev].element[ends*(node1-1)]].size[2];
           x2=E->ECO[lev][E->NEI[lev].element[(node2-1)*ends]].size[2];

           E->Interp[lev][2][node1] = x1;
           E->Interp[lev][2][node2] = x2;
           }

  if (3==dims)
        for(i=1;i<=nox;i++)   /* Fill in gaps in y direction */
            for(j=1;j<=noz;j++)
                for(k=2;k<noy;k+=2)   {
                    node0 = j + (i-1)*noz + (k-1)*nox*noz; /* this node */
                    node1 = node0 - nox*noz;
                    node2 = node0 + nox*noz;
                    x1=E->ECO[lev][E->NEI[lev].element[ends*(node1-1)]].size[3];
                    x2=E->ECO[lev][E->NEI[lev].element[(node2-1)*ends]].size[3];


                    E->Interp[lev][3][node1] = x1;
                    E->Interp[lev][3][node2] = x2;
                    }
return;
}


/* ==================================================== */
void inject_scalar_e(E,start_lev,AU,AD)

     struct All_variables *E;
     int start_lev;
     float *AU,*AD;  /* data on upper/lower mesh  */
{
    int i,j,m;
    int el,node,e;
    float average,w;

    const int sl_minus = start_lev-1;
    const int nels_minus=E->mesh.NEL[start_lev-1];
    const int dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
    const int n_minus=nels_minus*vpts;

    for(i=1;i<=n_minus;i++)
       AD[i] = 0.0;

    for(el=1;el<=nels_minus;el++)
        for(i=1;i<=ends;i++)                {
            e = E->EL[sl_minus][el].sub[i];
            AD[(el-1)*vpts+i] = AU[e];
            }

return;
}

void project_scalar(E,start_lev,AU,AD)

     struct All_variables *E;
     int start_lev;
     float *AU,*AD;  /* data on upper/lower mesh  */
{
    int i,j,m;
    int el,node,node1;
    float average,w;

    const int sl_minus = start_lev-1;
    const int nno_minus=E->mesh.NNO[start_lev-1];
    const int nels_minus=E->mesh.NEL[start_lev-1];
    const int  dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];
    const double weight=(double) 1.0/ends;

   for(i=1;i<=nno_minus;i++)
     AD[i] = 0.0;

        for(el=1;el<=nels_minus;el++)
            for(i=1;i<=ends;i++) {
                average=0.0;
                node1 = E->EL[sl_minus][el].sub[i];
                for(j=1;j<=ends;j++)                     {
                    node=E->IEN[start_lev][node1].node[j];
                    average += AU[node];
                    }

                w=weight*average;

                node= E->IEN[sl_minus][el].node[i];

                AD[node] += w * E->TWW[sl_minus][el].node[i];
         }

     for(i=1;i<=nno_minus;i++)  {
       AD[i] *= E->MASS[sl_minus][i];
       }


return;
}



/* ==================================================== */
void project_scalar_e(E,start_lev,AU,AD)

     struct All_variables *E;
     int start_lev;
     float *AU,*AD;  /* data on upper/lower mesh  */
{
    int i,j,m;
    int el,node,e;
    float average,w;

    const int sl_minus = start_lev-1;
    const int nels_minus=E->mesh.NEL[start_lev-1];
    const int  dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];
    const double weight=(double) 1.0/ends;
    const int vpts=vpoints[E->mesh.nsd];
    const int n_minus=nels_minus*vpts;

    for(i=1;i<=n_minus;i++)
       AD[i] = 0.0;

        for(el=1;el<=nels_minus;el++)    {
            average=0.0;
            for(i=1;i<=ends;i++) {
                e = E->EL[sl_minus][el].sub[i];
                average += log10((double)AU[e]);
                }
             
            AD[el] = pow((double)(10),(average*weight)); 
            }
return;
}

