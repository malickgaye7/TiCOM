/* Functions to assemble the element k matrices and the element f vector.
   Note that for the regular grid case the calculation of k becomes repetitive 
   to the point of redundancy. */

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <sys/time.h>
#include <sys/resource.h>


static double Dl[5][5] =   
	{	{ 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 1.0, 1.0, 0.0, 1.0 },
		{ 0.0, 1.0, 1.0, 0.0, 1.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 1.0, 1.0, 0.0, 1.0 } };

static double Dm[5][5] =   
	{	{ 0.0, 0.0, 0.0, 0.0, 0.0 },
		{ 0.0, 2.0, 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 2.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 1.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 2.0 } };



/* ================================================================
   Function to assemble the global  F vector.
                     +
   Function to get the global H vector (mixed method driving terms)
   ================================================================ */

void assemble_forces(E,penalty)
     struct All_variables *E;
     int penalty;
{
  double elt_f[24],elt_h[1]; 
  int el,p,i,a,a1,a2,a3,e,ii,jj,kk,elx,ely,elz,node,temp_dims;
  FILE *fp;
  char output_file[255];

  void get_elt_f();
  void matrix_transform_force();
  void get_elt_h();
  void strip_bcs_from_residual();
 
  double vdot();

  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[E->mesh.nsd];
  const int neq=E->mesh.neq;
  const int npno=E->mesh.npno;
  const int nel=E->mesh.nel;
  const int lev=E->mesh.levmax;

  for(a=0;a<neq;a++)
    E->F[a] = 0.0;
  
  /*for(a=0;a<npno;a++)
    E->H[a]=0.0;*/

   for (e=1;e<=nel;e++)  {

     get_elt_f(E,e,elt_f,penalty,1);
	  /*get_elt_h(E,e,elt_h,penalty);
	    E->H[e] = elt_h[0];  /* due to single pressure node per element */

     for(a=1;a<=ends;a++)          {
       a1=E->lm[e].node[a].doff[1];
       p=(a-1)*dims;
       E->F[a1] += elt_f[p];
       a2=E->lm[e].node[a].doff[2];
       E->F[a2] += elt_f[p+1];
       if (dims==3)   {
         a3=E->lm[e].node[a].doff[3];
         E->F[a3] += elt_f[p+2];
         }  
       }  

     }

   strip_bcs_from_residual(E,E->F,lev);

  return;
  }



/*==============================================================
  Function to supply the element k matrix for a given element e.
  ==============================================================  */

void get_elt_k(E,el,elt_k,level,penalty)
     struct All_variables *E;
     int el;
     double elt_k[24*24];
     int level;
     int penalty;
{
    double bdbmu[4][4];
    double bdbl[4][4];
    
    void print_elt_k();
    
    int pn,qn,ad,bd,k;

    int a,b,i,j,p,q,nint;
    double W[9],RM2[9],RMP[9],r;
    void get_global_shape_fn();
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    double xk[3][5],visc[9],temp;

 
    double Ba[5][3],Bb[5][3]; /* used for axi only */

    const int n=loc_mat_size[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
    const int ppts=ppoints[E->mesh.nsd];
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;

       
    get_global_shape_fn(E,el,&GN,&GNx,&dOmega,xk,penalty,level); 
				/* no pressure terms if not penalty */

    /* Note N[a].gauss_pt[n] is the value of shape fn a at the nth gaussian
       quadrature point. Nx[d] is the derivative wrt x[d]. */
       
    for(nint=1;nint<=vpts;nint++) {
	W[nint] = g_point[nint].weight[dims-1] * dOmega.vpt[nint] * E->EVI[level][(el-1)*vpts+nint];
    }  
 

     for(a=1;a<=ends;a++)
	    for(b=a;b<=ends;b++)   {
		ad=dims*(a-1);
		bd=dims*(b-1);
		bdbmu[1][1]=bdbmu[1][2]=
	        bdbmu[2][1]=bdbmu[2][2]=0.0;
	    
		if(dims==2)
		    for(k=1;k<=VPOINTS2D;k++)     {
			bdbmu[1][1] += W[k]*(
                             2.0*GNx.vpt[GNVXINDEX(0,a,k)]*GNx.vpt[GNVXINDEX(0,b,k)]
                            +(GNx.vpt[GNVXINDEX(1,a,k)]-xk[2][k]*E->N.vpt[GNVINDEX(a,k)])*
                             (GNx.vpt[GNVXINDEX(1,b,k)]-xk[2][k]*E->N.vpt[GNVINDEX(b,k)]) 
                            +2.0*xk[1][k]*xk[1][k]*xk[2][k]*xk[2][k]*
                             E->N.vpt[GNVINDEX(a,k)]*E->N.vpt[GNVINDEX(b,k)] );

			bdbmu[1][2] += W[k]*(
                             2.0*xk[2][k]*GNx.vpt[GNVXINDEX(0,a,k)]*E->N.vpt[GNVINDEX(b,k)]
                            +(GNx.vpt[GNVXINDEX(1,a,k)]-xk[2][k]*E->N.vpt[GNVINDEX(a,k)])*
                              GNx.vpt[GNVXINDEX(0,b,k)]+2.0*xk[1][k]*xk[2][k]*xk[2][k]* 
                              E->N.vpt[GNVINDEX(a,k)]*E->N.vpt[GNVINDEX(b,k)] );

			bdbmu[2][1] += W[k]*(
                             2.0*xk[2][k]*E->N.vpt[GNVINDEX(a,k)]*GNx.vpt[GNVXINDEX(0,b,k)]
                            +GNx.vpt[GNVXINDEX(0,a,k)]*
                            (GNx.vpt[GNVXINDEX(1,b,k)]-xk[2][k]*E->N.vpt[GNVINDEX(b,k)]) 
                            +2.0*xk[1][k]*xk[2][k]*xk[2][k]*E->N.vpt[GNVINDEX(a,k)]*E->N.vpt[GNVINDEX(b,k)] );

			bdbmu[2][2] += W[k]*(
                             2.0*xk[2][k]*xk[2][k]*E->N.vpt[GNVINDEX(a,k)]*E->N.vpt[GNVINDEX(b,k)]
                            +2.0*GNx.vpt[GNVXINDEX(1,a,k)]*GNx.vpt[GNVXINDEX(1,b,k)]
                            +GNx.vpt[GNVXINDEX(0,a,k)]*GNx.vpt[GNVXINDEX(0,b,k)] 
                            +2.0*xk[2][k]*xk[2][k]*E->N.vpt[GNVINDEX(a,k)]*E->N.vpt[GNVINDEX(b,k)] );
		    }
	    
		pn=ad*n+bd;
		qn=bd*n+ad;
	    
		elt_k[pn] = bdbmu[1][1] ; /* above */
		elt_k[pn+1] = bdbmu[1][2] ;	
		elt_k[pn+n] = bdbmu[2][1] ;
		elt_k[pn+n+1] = bdbmu[2][2] ;
	
		elt_k[qn] = bdbmu[1][1] ; /* below diag */
		elt_k[qn+n] = bdbmu[1][2] ;
		elt_k[qn+1] = bdbmu[2][1] ;
		elt_k[qn+n+1] = bdbmu[2][2] ;
	
	    } /*  Sum over all the a,b's to obtain full  elt_k matrix */
     

    return; 
}


/* =============================================
   General calling function for del_squared: 
   according to whether it should be element by
   element or node by node.
   ============================================= */

void assemble_del2_u(E,u,Au,level,strip_bcs)
     struct All_variables *E;
     double *u,*Au;
     int level;
     int strip_bcs;
{
  void e_assemble_del2_u();
  void n_assemble_del2_u();
 
  if(E->control.NMULTIGRID||E->control.NASSEMBLE)
    n_assemble_del2_u(E,u,Au,level,strip_bcs);
  else
    e_assemble_del2_u(E,u,Au,level,strip_bcs);
 
  return;
}

/* ======================================
   Assemble del_squared_u vector el by el
   ======================================   */

void e_assemble_del2_u(E,u,Au,level,strip_bcs)
  struct All_variables *E;
  double *u,*Au;
  int level;
  int strip_bcs;

{ 
  int  el,e,i,a,b,a1,a2,a3,ii;
  double elt_k[24*24],U[24],AU[24],alpha=1.0,beta=0.0;
  
  int indx[24];
  char uplo='U';
    
  void strip_bcs_from_residual();
  void get_elt_k();
   
  double U1[24];
   
  const int n=loc_mat_size[E->mesh.nsd];
  const int ends=enodes[E->mesh.nsd]; 
  const int dims=E->mesh.nsd,dofs=E->mesh.dof; 
  const int nel=E->mesh.NEL[level];
  const int neq=E->mesh.NEQ[level];

  for(i=0;i<neq;i++)
     Au[i] = 0.0;
 
  for(e=1;e<=nel;e++)   {

     for(a=1;a<=ends;a++) {
	a1 = E->LMD[level][e].node[a].doff[1];   
	a2 = E->LMD[level][e].node[a].doff[2];   
	if(dims ==3 ) a3 = E->LMD[level][e].node[a].doff[3];   
	for(b=1;b<=ends;b++) {
	    ii = (a*n+b)*dims-(dims*n+dims);
	    if(dims==2)  {
			/* i=1, j=1,2 */ 
		Au[a1] += 
	                E->elt_k[level][e].k[ii  ] * 
			u[E->LMD[level][e].node[b].doff[1]] 
	              + E->elt_k[level][e].k[ii+1] *
		        u[E->LMD[level][e].node[b].doff[2]]; 
			/* i=2, j=1,2 */
		Au[a2] += 
		        E->elt_k[level][e].k[ii+n  ] * 
			u[E->LMD[level][e].node[b].doff[1]] 
		      + E->elt_k[level][e].k[ii+n+1] * 
			u[E->LMD[level][e].node[b].doff[2]]; 
       	        }
	    if(3==dims)   {
	          /* i=1, j=1,2,3 */ 
		Au[a1] += 
		        E->elt_k[level][e].k[ii] * 
			u[E->LMD[level][e].node[b].doff[1]]   
		      + E->elt_k[level][e].k[ii+1] * 
			u[E->LMD[level][e].node[b].doff[2]]   
		      + E->elt_k[level][e].k[ii+2] * 
			u[E->LMD[level][e].node[b].doff[3]]; 
		/* i=2, j=1,2,3 */
		Au[a2] += 
		        E->elt_k[level][e].k[ii+n] * 
			u[E->LMD[level][e].node[b].doff[1]] 
		      + E->elt_k[level][e].k[ii+n+1] * 
			u[E->LMD[level][e].node[b].doff[2]]   
		      + E->elt_k[level][e].k[ii+n+2] * 
			u[E->LMD[level][e].node[b].doff[3]]; 
		/* i=3, j=1,2,3 */
		Au[a3] += 
		        E->elt_k[level][e].k[ii+n+n] * 
			u[E->LMD[level][e].node[b].doff[1]] 
		      + E->elt_k[level][e].k[ii+n+n+1] * 
			u[E->LMD[level][e].node[b].doff[2]]   
		      + E->elt_k[level][e].k[ii+n+n+2] * 
			u[E->LMD[level][e].node[b].doff[3]]; 
	   
	        }    /* end for if dims==3 */
 	    }         /* end for loop b */
        }             /* end for loop a */

     }

  if(strip_bcs)
     strip_bcs_from_residual(E,Au,level);

  return; }


/* ======================================================
   Assemble Au using stored, nodal coefficients.
   ====================================================== */

void n_assemble_del2_u(E,u,Au,level,strip_bcs)
     struct All_variables *E;
     double *u,*Au;
     int level;
     int strip_bcs;
{
    int node, e,i,max_eqn;
    int eqn1,eqn2,eqn3,loc0,loc1,loc2,loc3;
   
    double U1,U2,U3;
    void strip_bcs_from_residual();

    static int been_here=0;
 
    int *C;
    higher_precision *B1,*B2,*B3;

    const int neq=E->mesh.NEQ[level];
    const int nno=E->mesh.NNO[level];
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  
    for(e=0;e<=neq;e++)
	Au[e]=0.0;

    loc0 = 1;

    if(3==dims)
	for(e=1;e<=nno;e++)     {

	    if(E->NODE[level][e] & OFFSIDE) 
		continue;
	     eqn1=E->ID[level][e].doff[1];
	     eqn2=E->ID[level][e].doff[2];
	     eqn3=E->ID[level][e].doff[3];

  	     U1 = u[eqn1];
	     U2 = u[eqn2];
	     U3 = u[eqn3];

	     max_eqn = E->Node_eqn[level][e+1]-E->Node_eqn[level][e];
 
	     C=E->Node_map[level] + E->Node_eqn[level][e];
	     B1=E->Eqn_k[level]+E->Node_k_id[level][e];
	     B2=E->Eqn_k[level]+E->Node_k_id[level][e]+max_eqn;
	     B3=E->Eqn_k[level]+E->Node_k_id[level][e]+2*max_eqn;
	   
	     for(i=0;i<max_eqn;i++) {
	 	  Au[C[i]] += B1[i]*U1+B2[i]*U2+B3[i]*U3;
	          }

	     }
    else
	for(e=1;e<=nno;e++)     {

	    if(E->NODE[level][e] & OFFSIDE) 
		continue;
	    eqn1=E->ID[level][e].doff[1];
	    eqn2=E->ID[level][e].doff[2];
	  
	    U1 = u[eqn1];
	    U2 = u[eqn2]; 
	    
	    max_eqn = E->Node_eqn[level][e+1]-E->Node_eqn[level][e];

	    C=E->Node_map[level] + E->Node_eqn[level][e];
	    B1=E->Eqn_k[level]+E->Node_k_id[level][e];
	    B2=E->Eqn_k[level]+E->Node_k_id[level][e]+max_eqn;


	    for(i=0;i<max_eqn;i++) {
     	        Au[C[i]] += B1[i]*U1+B2[i]*U2;
		}

	    }

    if (strip_bcs)
	strip_bcs_from_residual(E,Au,level);

    been_here++;

    return; 
}


void build_diagonal_of_K(E,el,elt_k,level)
     struct All_variables *E;
     int level,el;
     double elt_k[24*24];

{
    int a,a1,a2,p;
     
    const int n=loc_mat_size[E->mesh.nsd];
    const int dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];

    for(a=1;a<=ends;a++) {
	    /* dirn 1 */
	    a1 = E->LMD[level][el].node[a].doff[1];
	    p=(a-1)*dims; 
	    E->BI[level][a1] += elt_k[p*n+p];  

	    /* dirn 2 */
	    a2 = E->LMD[level][el].node[a].doff[2];
	    p=(a-1)*dims+1; 
	    E->BI[level][a2] += elt_k[p*n+p];  

	    /* dirn 3 */
	    if(3==dims){
	        a1 = E->LMD[level][el].node[a].doff[3];
		p=(a-1)*dims+2; 
	        E->BI[level][a1] += elt_k[p*n+p];  
	        }
            } 	    

  return;
}
  
void build_diagonal_of_Ahat(E,level)
    struct All_variables *E;
    int level;
{
    double assemble_dAhatp_entry();
 
    double BU;
    int e,npno,neq;
    float time,time0,CPU_time();
   
    npno = E->mesh.NPNO[level];
    neq=E->mesh.NEQ[level];
       
    for(e=1;e<=npno;e++)
	E->BPI[level][e]=1.0;

    if(!E->control.precondition)
	return;
    
    for(e=1;e<=npno;e++)  {
	BU=assemble_dAhatp_entry(E,e,level);
	if(BU != 0.0)
	    E->BPI[level][e] = 1.0/BU;
	else
	    E->BPI[level][e] = 1.0;
  }

    return;
}

/* ==========================================
   Assemble a div_u vector element by element
   ==========================================  */

void assemble_div_u(E,U,divU,level)
     struct All_variables *E;
     double *U,*divU;
     int level;       
{ 
    int e,j1,j2,j3,p,a,b;
    higher_precision elt_g[24][1];
    void get_elt_g();
    
    const int nel=E->mesh.NEL[level];
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int npno=E->mesh.NPNO[level];
    
    for(e=1;e<=npno;e++)
	divU[e] = 0.0;
   
    if (dims==3) 
       for(a=1;a<=ends;a++)   {
	  p = (a-1)*dims;
          for(e=1;e<=nel;e++) {
	    j1= E->LMD[level][e].node[a].doff[1];
	    j2= E->LMD[level][e].node[a].doff[2];
	    j3= E->LMD[level][e].node[a].doff[3];
	    /* for(b=0;b<ploc_mat_size[E->mesh.nsd];b++) */
	    divU[e] += E->elt_del[level][e].g[p  ][0] * U[j1]
	             + E->elt_del[level][e].g[p+1][0] * U[j2]
	             + E->elt_del[level][e].g[p+2][0] * U[j3];
	    }
	 } 
    else if (dims==2) 
       for(a=1;a<=ends;a++)   {
	  p = (a-1)*dims;
          for(e=1;e<=nel;e++) {
	    j1= E->LMD[level][e].node[a].doff[1];
	    j2= E->LMD[level][e].node[a].doff[2];
	    /* for(b=0;b<ploc_mat_size[E->mesh.nsd];b++) */
	    divU[e] += E->elt_del[level][e].g[p  ][0] * U[j1]
	             + E->elt_del[level][e].g[p+1][0] * U[j2];
	    }
	 } 
    
    return;
}


/* ==========================================
   Assemble a grad_P vector element by element
   ==========================================  */

void assemble_grad_p(E,P,gradP,lev)
     struct All_variables *E;
     double *P,*gradP;
     int lev;
    
{
  int el,e,i,j1,j2,p,a,nel,neq;
  void strip_bcs_from_residual();
  higher_precision elt_g[24][1];
  void get_elt_g();
 
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd,dofs=E->mesh.dof;

  nel=E->mesh.NEL[lev];
  neq=E->mesh.NEQ[lev];

  for(i=0;i<=neq;i++)
    gradP[i] = 0.0;
 
  for(e=1;e<=nel;e++) {
     
	if(0.0==P[e])
	    continue;

	for(a=1;a<=ends;a++)       {
	     p = (a-1)*dims;
             j1=E->LMD[lev][e].node[a].doff[1];
             j2=E->LMD[lev][e].node[a].doff[2];
		        /*for(b=0;b<ploc_mat_size[E->mesh.nsd];b++)  */ 
             gradP[j1] += E->elt_del[lev][e].g[p][0] * P[e];
             gradP[j2] += E->elt_del[lev][e].g[p+1][0] * P[e];
	     if(dims==3) {
               j1=E->LMD[lev][e].node[a].doff[3];
               gradP[j1] += E->elt_del[lev][e].g[p+2][0] * P[e];
	       }
	     }

        }

  strip_bcs_from_residual(E,gradP,lev);

return; 
}

double assemble_dAhatp_entry(E,e,level)
     struct All_variables *E;
     int e,level;
    
{ 
    int i,j,p,a,b,node,ee,element,lnode,npno;
    void strip_bcs_from_residual();
    higher_precision elt_g[24][1];
    void get_elt_g();

    double gradP[81],divU;

    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;

    npno=E->mesh.NPNO[level];
    
    for(i=0;i<81;i++)
	gradP[i] = 0.0;
  
    divU=0.0;
 
    for(a=1;a<=ends;a++) {
      p = (a-1)*dims;
      j=E->LMD[level][e].node[a].doff[1];
      gradP[p] += E->BI[level][j]*E->elt_del[level][e].g[p][0];

      j=E->LMD[level][e].node[a].doff[2];
      gradP[p+1] += E->BI[level][j]*E->elt_del[level][e].g[p+1][0];
	    
      if(3==dims) { 
        j=E->LMD[level][e].node[a].doff[3];
        gradP[p+2] += E->BI[level][j]*E->elt_del[level][e].g[p+2][0];
	    }
	  }

   
    /* calculate div U from the same thing .... */

    /* only need to run over nodes with non-zero grad P, i.e. the ones in
       the element accessed above, BUT it is only necessary to update the
       value in the original element, because the diagonal is all we use at
       the end ... */

    for(b=1;b<=ends;b++) {
      p = (b-1)*dims;	   
      divU +=E->elt_del[level][e].g[p][0] * gradP[p];	    
      divU +=E->elt_del[level][e].g[p+1][0] * gradP[p+1];
      if(3==dims)
        divU +=E->elt_del[level][e].g[p+2][0] * gradP[p+2];
      }
   
return(divU);  }


/*==============================================================
  Function to supply the element g matrix for a given element e.
  ==============================================================  */

void get_elt_g(E,el,elt_del,level)
     struct All_variables *E;
     int el;
     higher_precision elt_del[24][1];
     int level;

{  void get_global_shape_fn();
   int p,a,nint;
   double dGNdash[3];
   double aaa,xk[3][5];
   int lmsize;

   struct Shape_function GN;
   struct Shape_function_dA dOmega;
   struct Shape_function_dx GNx;
    
   const int dims=E->mesh.nsd,dofs=E->mesh.dof;
   const int ends=enodes[dims];
   const int vpts=vpoints[dims];
  
   get_global_shape_fn(E,el,&GN,&GNx,&dOmega,xk,2,level);

   aaa=p_point[1].weight[dims-1] * dOmega.ppt[1];
  
   for(a=1;a<=ends;a++)  {
       p=dims*(a-1);
       elt_del[p  ][0]=-(GNx.ppt[GNPXINDEX(0,a,1)]+xk[2][1]*xk[1][1]*E->N.ppt[GNPINDEX(a,1)])*aaa;
       elt_del[p+1][0]=-(GNx.ppt[GNPXINDEX(1,a,1)]+xk[2][1]*2.0*E->N.ppt[GNPINDEX(a,1)])*aaa;
     }
       
   return; 
 }




/* ===============================================================
   Function to create the element pressure-forcing vector (due
   to imposed velocity boundary conditions, mixed method).
   =============================================================== */

void get_elt_h(E,el,elt_h,penalty)
     struct All_variables *E;
     int el;
     double elt_h[1];
     int penalty;
{
    int aid,i,p,a,b,d,j,k,q,global,got_g;
    unsigned int type;
    double elt_g[24][1];
    void get_elt_g();

    for(p=0;p<1;p++) elt_h[p] = 0.0;
   
    if (penalty) return; /* no h term at all */

    got_g = 0;
  
  type=VBX;
  for(i=1;i<=E->mesh.nsd;i++)
    { for(a=1;a<=enodes[E->mesh.nsd];a++)
	{ if (E->node[E->ien[el].node[a]] & type)
	    { if(!got_g) 
		{  get_elt_g(E,el,elt_g,E->mesh.levmax);
		   got_g++;
		 }
	      
	      p=E->mesh.nsd*(a-1) + i - 1;
	      for(b=1;b<=pnodes[E->mesh.nsd];b++)
		{ q = b-1; 
		  elt_h[q] -= elt_g[p][q] * E->VB[i][E->ien[el].node[a]];
		}
	    }
	}
      type *= (unsigned int) 2;
    }
   return;
}
	
/*=================================================================
  Function to create the element force vector (allowing for b.c.'s)
  ================================================================= */

void get_elt_f(E,el,elt_f,penalty,bcs)
     struct All_variables *E;
     int el;
     double elt_f[24];
     int penalty,bcs;
   
{
  /* 	Three parts to the element force vector:
	. Globally derived imposed forcing terms.
	. Non-zero velocity boundary conditions (propogated to interior of element)
	. Locally applicable stress boundary conditions.
	This follows the treatment given in Hughes.                   */

  int aid,i,p,a,b,d,j,k,q;
  int node[5],back_front,got_elt_k,nodea,nodeb;
  unsigned int type;

  double xk[3][5],force[9],force_at_gs[9],stress[9],elt_k[24*24];
  double vector[4],magnitude;
  double tmp;
  
  void get_global_shape_fn();
  void get_global_1d_shape_fn();

  struct Shape_function GN;
  struct Shape_function_dA dOmega;
  struct Shape_function_dx GNx;
  struct Shape_function1 GM;
  struct Shape_function1_dA dGammax;
  
  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int n=loc_mat_size[dims];
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];

  get_global_shape_fn(E,el,&GN,&GNx,&dOmega,xk,0,E->mesh.levmax);
  get_global_1d_shape_fn(E,el,&GM,&dGammax,0);
  

  for(p=0;p<n;p++) elt_f[p] = 0.0;
  
  for(i=1;i<=dims;i++)  {
    if(i!=2) 
	  for(p=1;p<=ends;p++)
	    force[p] = 0.0;
    else
	  for(p=1;p<=ends;p++)
	    force[p] = E->buoyancy[E->ien[el].node[p]];
      
    for(a=1;a<=ends;a++)  {
	  nodea=E->ien[el].node[a];
	  p= dims*(a-1)+i-1;
	  
	  for(j=1;j<=vpts;j++) {   /*compute force at each int point */
		 force_at_gs[j] = 0.0;
		 for(k=1;k<=ends;k++)
		     force_at_gs[j] += force[k] * E->N.vpt[GNVINDEX(k,j)] ;
		 }

	  for(j=1;j<=vpts;j++)     /*compute sum(Na(j)*F(j)*det(j)) */
		 elt_f[p] += force_at_gs[j] * E->N.vpt[GNVINDEX(a,j)]
			      *dOmega.vpt[j]*g_point[j].weight[dims-1] ;
		   
	  /* imposed velocity terms */

	  if(bcs)  {
	     got_elt_k = 0; 
	     type=VBX;
	     for(j=1;j<=dims;j++) {
		   for(b=1;b<=ends;b++) {
		      nodeb=E->ien[el].node[b];
		      if ((E->node[nodeb] & type) && (E->VB[j][nodeb] != 0.0)){
		 	    if(!got_elt_k) {
			      get_elt_k(E,el,elt_k,E->mesh.levmax,penalty);
			      got_elt_k = 1;
			      }
			    q = dims*(b-1)+j-1;
			    if(p!=q) {
			      elt_f[p] -= elt_k[p*n+q] * E->VB[j][nodeb];
/* printf("el %d: dirn=%d, vbc found at node %d of %g\n",el,j,b,E->VB[j][nodeb]);*/
			      }
		        }
		      }  /* end for b */
		   type *= (unsigned int) 2;
	       }      /* end for j */
	     }      /* end if for if bcs */
		 
	   }
    } /*  Complete the loops for a,i  	*/



  return; 
}

/* =================================================================
 subroutine to get augmented lagrange part of stiffness matrix
================================================================== */

void get_aug_k(E,el,elt_k,level,penalty)
     struct All_variables *E;
     int el;
     double elt_k[24*24];
     int level;
     int penalty;
{
     int i,j,k,p[9],a,b,nodea,nodeb;
     double Visc;

     const int n=loc_mat_size[E->mesh.nsd];
     const int ends=enodes[E->mesh.nsd];
     const int vpts=vpoints[E->mesh.nsd];
     const int dims=E->mesh.nsd;

     Visc = 0.0;
     for(a=1;a<=vpts;a++) {
	  p[a] = (a-1)*dims;
	  Visc += E->EVI[level][(el-1)*vpts+a];
       }
     Visc = Visc/vpts;

     for(a=1;a<=ends;a++) {
        nodea=E->IEN[level][el].node[a];
        for(b=1;b<=ends;b++) {
           nodeb=E->IEN[level][el].node[b];      /* for Kab dims*dims  */
	   i = (a-1)*n*dims+(b-1)*dims;
	   elt_k[i  ] += Visc*E->control.augmented*
	              E->elt_del[level][el].g[p[a]][0]*
		      E->elt_del[level][el].g[p[b]][0];   /*for 11 */
	   elt_k[i+1] += Visc*E->control.augmented*
	              E->elt_del[level][el].g[p[a]][0]*
		      E->elt_del[level][el].g[p[b]+1][0];  /* for 12 */
	   elt_k[i+n] += Visc*E->control.augmented*          
	              E->elt_del[level][el].g[p[a]+1][0]*
		      E->elt_del[level][el].g[p[b]][0];    /* for 21 */
	   elt_k[i+n+1] += Visc*E->control.augmented*
	              E->elt_del[level][el].g[p[a]+1][0]*
		      E->elt_del[level][el].g[p[b]+1][0];  /* for 22 */

           if(3==dims) {
	       elt_k[i+2] += Visc*E->control.augmented*
	              E->elt_del[level][el].g[p[a]][0]*
		      E->elt_del[level][el].g[p[b]+2][0];  /* for 13 */
	       elt_k[i+n+2] += Visc*E->control.augmented*
	              E->elt_del[level][el].g[p[a]+1][0]*
		      E->elt_del[level][el].g[p[b]+2][0];  /* for 23 */
	       elt_k[i+n+n] += Visc*E->control.augmented*
	              E->elt_del[level][el].g[p[a]+2][0]*
		      E->elt_del[level][el].g[p[b]][0];    /* for 31 */
	       elt_k[i+n+n+1] += Visc*E->control.augmented*
	              E->elt_del[level][el].g[p[a]+2][0]*
		      E->elt_del[level][el].g[p[b]+1][0];  /* for 32 */
	       elt_k[i+n+n+2] += Visc*E->control.augmented*
	              E->elt_del[level][el].g[p[a]+2][0]*
		      E->elt_del[level][el].g[p[b]+2][0];  /* for 33 */
               }
           }
       }

   return;
   }



