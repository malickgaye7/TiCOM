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

/*========================================================
  Function to make the IEN array for a mesh of given 
  dimension. IEN is an externally defined structure array

  NOTE: this is not really general enough for new elements:
  it should be done through a pre-calculated lookup table.
  ======================================================== */

void construct_ien(E)
     struct All_variables *E;

{	
  int lev,p,q,r,rr,e1,e2,i,a,node,node2,e;
  int element,start,start1,nel,nno;
  int elz,elx,ely,nox,noy,noz,nozx,nozx1;


  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];

  for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
    elz = E->mesh.ELZ[lev];
    ely = E->mesh.ELY[lev];
    noz = E->mesh.NOZ[lev];
    noy = E->mesh.NOY[lev];
    nox = E->mesh.NOX[lev];
    elx = E->mesh.ELX[lev];
    nel=E->mesh.NEL[lev];
    nno=E->mesh.NNO[lev];

      for(r=1;r<=ely;r++) 
        for(q=1;q<=elx;q++)    
          for(p=1;p<=elz;p++)     {
		    element = (r-1)*elx*elz + (q-1)*elz  + p;	
		    start = (r-1)*noz*nox + (q-1)*noz + p;
		    for(rr=1;rr<=ends;rr++)
		         E->IEN[lev][element].node[rr]= start 
			     + offset[rr].vector[1]
			     + offset[rr].vector[0]*noz
			     + offset[rr].vector[2]*noz*nox;   
	        }    

    for(i=1;i<=nno;i++)
      E->NEI[lev].nels[i] = 0;

    for(e=1;e<=nel;e++)
      for(a=1;a<=ends;a++) {
        node=E->IEN[lev][e].node[a];
	    E->NEI[lev].nels[node]++;
        E->NEI[lev].element[(node-1)*ends+E->NEI[lev].nels[node]-1] = e;
	    E->NEI[lev].lnode[(node-1)*ends+E->NEI[lev].nels[node]-1] = a;
        }

     }     /* end loop for lev */

             /*  determine surface things */
    e = 0;
    for(element=1;element<=E->mesh.NEL[E->mesh.levmax];element++)
      if ( (element-1)%E->mesh.elz==0) {
        e ++;
        E->sien[e].node[1] = (E->ien[element].node[1]-1)/E->mesh.noz+1;
        E->sien[e].node[2] = (E->ien[element].node[4]-1)/E->mesh.noz+1;
        if (dims==3)   {
          E->sien[e].node[3] = (E->ien[element].node[8]-1)/E->mesh.noz+1;
          E->sien[e].node[4] = (E->ien[element].node[5]-1)/E->mesh.noz+1;
          }
        E->surf_element[e] = element;
        }
    E->mesh.snel = e;
    for (i=1;i<=E->mesh.nsf;i++)
      E->surf_node[i] = i*E->mesh.noz;


if (E->control.verbose) 
  for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
     fprintf(E->fp,"output_IEN_arrays %d\n",lev);
     if(dims==2) 
        for (i=1;i<=E->mesh.NEL[lev];i++)
           fprintf(E->fp,"%d %d %d %d %d\n",i,E->IEN[lev][i].node[1],E->IEN[lev][i].node[2],E->IEN[lev][i].node[3],E->IEN[lev][i].node[4]);
     else if(dims==3) 
        for (i=1;i<=E->mesh.NEL[lev];i++)
           fprintf(E->fp,"%d %d %d %d %d %d %d %d %d\n",i,E->IEN[lev][i].node[1],E->IEN[lev][i].node[2],E->IEN[lev][i].node[3],E->IEN[lev][i].node[4],E->IEN[lev][i].node[5],E->IEN[lev][i].node[6],E->IEN[lev][i].node[7],E->IEN[lev][i].node[8]);
     }    

  return;
}

/*============================================
  Function to make the ID array for above case
  ============================================ */

void construct_id(E)
     struct All_variables *E;
{ 
    int i,j,k,kk,i1,i2,j1,j2,k1,k2;	
    int eqn_count,node,eqn_countd;
    unsigned int type,doff;
    int lev, temp_dims;
    int nox,noy,noz;
    int elx,ely,elz,rr,element;

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];
    const int addi_dof=additional_dof[dims];

    for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
      eqn_count = 0;
      elz = E->mesh.ELZ[lev];
      ely = E->mesh.ELY[lev];
      noz=E->mesh.NOZ[lev];
      noy=E->mesh.NOY[lev];

      for(node=1;node<=E->mesh.NNO[lev];node++)
	  for(doff=1;doff<=dims;doff++)  {
              E->ID[lev][node].doff[doff] = eqn_count;
              eqn_count += 1;
              }
      E->mesh.NEQ[lev] = eqn_count;

      }


    E->mesh.neq = E->mesh.NEQ[E->mesh.levmax];  /*  Total NUMBER of independent variables  */


if (E->control.verbose)
  for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
      fprintf(E->fp,"output_ID_arrays %d\n",lev);
      if (dims==2)
        for (i=1;i<=E->mesh.NNO[lev];i++)
          fprintf(E->fp,"%d %d %d \n",i,E->ID[lev][i].doff[1],E->ID[lev][i].doff[2]);
      else if (dims==3)
        for (i=1;i<=E->mesh.NNO[lev];i++)
          fprintf(E->fp,"%d %d %d %d\n",i,E->ID[lev][i].doff[1],E->ID[lev][i].doff[2],E->ID[lev][i].doff[3]);
      }

    return; 
    }

/*==========================================================
  Function to construct  the LM array from the ID and IEN arrays 
  ========================================================== */

void construct_lm(E)
     struct All_variables *E;
{	
  int i,a,e;
  int lev,eqn_no;
  int nel, nel2;
  
  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];
  const int addi_dof=additional_dof[dims];
  

  for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)  {
    nel=E->mesh.NEL[lev];
    for(e=1;e<=nel;e++)
       for(a=1;a<=ends;a++)     {
	           E->LMD[lev][e].node[a].doff[1] = E->ID[lev][E->IEN[lev][e].node[a]].doff[1]; 
	           E->LMD[lev][e].node[a].doff[2] = E->ID[lev][E->IEN[lev][e].node[a]].doff[2]; 
	           if(dims==3)
                      E->LMD[lev][e].node[a].doff[3] = E->ID[lev][E->IEN[lev][e].node[a]].doff[3]; 
		       }
   }                /* end for level */
     
  if(E->control.verbose) 
     if(dims==3)
        for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)  {
          fprintf(E->fp,"output_LM_arrays %d\n",lev);
          nel=E->mesh.NEL[lev];
          for(e=1;e<=nel;e++)   {
            for(a=1;a<=ends;a++) 
	      fprintf(E->fp,"%d %d %d %d %d\n",e,a,E->LMD[lev][e].node[a].doff[1],E->LMD[lev][e].node[a].doff[2],E->LMD[lev][e].node[a].doff[3]); 
	      }
        }
     else if(E->mesh.nsd==2)
        for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)  {
          fprintf(E->fp,"output_LM_arrays %d\n",lev);
          nel=E->mesh.NEL[lev];
          for(e=1;e<=nel;e++)   
            for(a=1;a<=ends;a++) 
	      fprintf(E->fp,"%d %d %d %d\n",e,a,E->LMD[lev][e].node[a].doff[1],E->LMD[lev][e].node[a].doff[2]); 
        }
  fflush(E->fp);

 
  return;	
}


/* =====================================================
   Function to build the local node matrix indexing maps
   ===================================================== */

void construct_node_maps(E)
    struct All_variables *E;
{
    float initial_time,CPU_time();

    int el,n,nn,lev,i,j,jj;
    int node1,eqn1,loc1,count,found,element;
    int neq,nno,temp_dims;
   
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];
    const int max_eqn = max_eqn_interaction[dims];
    const int addi_dof = additional_dof[dims];

    for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)   {
       neq=E->mesh.NEQ[lev];
       nno=E->mesh.NNO[lev];
       E->Node_map[lev]=(int *) malloc ((nno+5)*max_eqn*sizeof(int));
        
       for(i=0;i<=(nno+1)*max_eqn;i++) 
	   E->Node_map[lev][i] = neq+1;  /* DANGER !!! */
       for(i=1;i<=(nno+1);i++) 
	   E->Node_eqn[lev][i] = 0;
      
       loc1 = 1;
       count = 0;
       for(nn=1;nn<=nno;nn++) {
	   if(E->NODE[lev][nn] & OFFSIDE)
	       continue;
	   for(el=1;el<=E->NEI[lev].nels[nn];el++)  {
	       element = E->NEI[lev].element[(nn-1)*ends+el-1]; 
	       for(n=1;n<=ends;n++) {
		   node1=E->IEN[lev][element].node[n]; /*global node number*/
		   for(i=1;i<=dims;i++) {
		       eqn1=E->ID[lev][node1].doff[i];
 
		       found=0;
		       for(j=0;j<=count;j++)
			   if(E->Node_map[lev][loc1+j] == eqn1) { /* found, index next equation */
			       found++;
			       break;
	             	       }
		       
		       if(! found) {
			   E->Node_map[lev][loc1+count] = eqn1;
			   count++;
   		           }
		       }      /* end i */
	           }          /* end n */
	       }              /* end el */
/* fprintf(E->fp,"%d %d %d \n",nn,E->NEI[lev].nels[nn],count); 
*/	   
           E->Node_eqn[lev][nn] = loc1;
	   loc1 += count;
	   count=0;
           }                      /* end for node */

       E->Node_eqn[lev][nno+1] = loc1;       /* mark the end of Node_eqn */

       loc1 = 0;
       for(nn=1;nn<=nno;nn++) {
	   if(E->NODE[lev][nn] & OFFSIDE)
	       continue;

	   E->Node_k_id[lev][nn] = loc1;
 	   loc1 += (E->Node_eqn[lev][nn+1]-E->Node_eqn[lev][nn]) * dims;
	   }

       E->Eqn_k[lev] = (higher_precision *) malloc ((loc1+5) * sizeof(higher_precision));

       E->mesh.matrix_size[lev] = loc1 + 1;
       }                          /* end for level */
    
    return;
}


void construct_node_ks(E)
     struct All_variables *E;
{ 
    int level,i,j,k,e;
    int node,node1,eqn1,eqn2,eqn3,loc0,loc1,loc2,loc3,found,element,index,pp,qq;
    int neq,nno,nel,max_eqn;
   
    double elt_K[24*24];
    double w1,w2,w3,ww1,ww2,ww3;
    higher_precision *B1,*B2,*B3;
   
    void get_elt_k();
    void matrix_transform_K();
    void get_aug_k();
    void build_diagonal_of_K();
   
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];
    const int lms=loc_mat_size[E->mesh.nsd];
 
    for(level=E->mesh.levmax;level>=E->mesh.levmin;level--)   {
	neq=E->mesh.NEQ[level];
	nel=E->mesh.NEL[level];
	nno=E->mesh.NNO[level];

	for(i=0;i<=(neq+1);i++) 
	    E->BI[level][i] = 0.0; 
    for(i=0;i<=E->mesh.matrix_size[level];i++) 
        E->Eqn_k[level][i] = 0.0; 
	
    for(element=1;element<=nel;element++) {

	    get_elt_k(E,element,elt_K,level,0);

	    if (E->control.augmented_Lagr)
	         get_aug_k(E,element,elt_K,level,0);

        build_diagonal_of_K(E,element,elt_K,level);
	     
	    for(i=1;i<=ends;i++) {  /* i, is the node we are storing to */
		node=E->IEN[level][element].node[i];
		if(E->NODE[level][node] & OFFSIDE)
		    continue;
		
		pp=(i-1)*dims;
		w1=w2=w3=1.0;
	
		loc0=E->Node_eqn[level][node];
		max_eqn = E->Node_eqn[level][node+1] - E->Node_eqn[level][node];
		loc1=E->Node_k_id[level][node];
		loc2=loc1+(E->LMD[level][element].node[i].doff[2]
			  -E->LMD[level][element].node[i].doff[1])*max_eqn;
		if (3==dims) loc3=loc1+(E->LMD[level][element].node[i].doff[3]
			         -E->LMD[level][element].node[i].doff[1])*max_eqn;
	
		if(E->NODE[level][node] & VBX) w1=0.0;
		if(E->NODE[level][node] & VBZ) w2=0.0;
		if(E->NODE[level][node] & VBY) w3=0.0;
		   
		for(j=1;j<=ends;j++) { /* j is the node we are receiving from */
		    ww1=ww2=ww3=1.0;
		    qq=(j-1)*dims;
		    node1=E->IEN[level][element].node[j];
		    eqn1=E->LMD[level][element].node[j].doff[1];
		    eqn2=E->LMD[level][element].node[j].doff[2];
		    if(3==dims) eqn3=E->LMD[level][element].node[j].doff[3];

		    if(E->NODE[level][node1] & VBX) ww1=0.0;
		    if(E->NODE[level][node1] & VBZ) ww2=0.0;
		    if(E->NODE[level][node1] & VBY) ww3=0.0; 
		    
		    /* search for direction 1*/

		    found=0;
		    for(k=0;k<max_eqn;k++)
			if(E->Node_map[level][loc0+k] == eqn1) { /* found, index next equation */
			    index=k;
			    found++;
			    break;
			}
		       
		    assert(found /* direction 1 */);

		    E->Eqn_k[level][loc1+index] +=  w1*ww1*elt_K[pp*lms+qq]; /* direction 1 */
		    E->Eqn_k[level][loc2+index] +=  w2*ww1*elt_K[(pp+1)*lms+qq]; /* direction 1 */
		    if(3==dims) E->Eqn_k[level][loc3+index] +=  w3*ww1*elt_K[(pp+2)*lms+qq]; /* direction 1 */
		    
		     /* search for direction 2*/

		    found=0;
		    for(k=0;k<max_eqn;k++)
			if(E->Node_map[level][loc0+k] == eqn2) { /* found, index next equation */
			    index=k;
			    found++;
			    break;
			}

		    assert(found /* direction 2 */);

		    E->Eqn_k[level][loc1+index] += w1*ww2*elt_K[pp*lms+qq+1]; /* direction 1 */
		    E->Eqn_k[level][loc2+index] += w2*ww2*elt_K[(pp+1)*lms+qq+1]; /* direction 2 */
		    if(3==dims) E->Eqn_k[level][loc3+index] += w3*ww2*elt_K[(pp+2)*lms+qq+1]; /* direction 3 */
		
		    /* search for direction 3*/
		   
		    if(3==dims) {
			found=0;
			for(k=0;k<max_eqn;k++)
			    if(E->Node_map[level][loc0+k] == eqn3) { /* found, index next equation */
				index=k;
				found++;
				break;
			    }
			

			assert(found /* direction 3 */);
			
			E->Eqn_k[level][loc1+index] += w1*ww3*elt_K[pp*lms+qq+2]; /* direction 1 */
			E->Eqn_k[level][loc2+index] += w2*ww3*elt_K[(pp+1)*lms+qq+2]; /* direction 2 */
			E->Eqn_k[level][loc3+index] += w3*ww3*elt_K[(pp+2)*lms+qq+2]; /* direction 3 */

		    }
		}
	    }
	}


        for(j=0;j<neq;j++) {
            if(E->BI[level][j] ==0.0)  fprintf(stderr,"level %d, equation %d/%d has zero diagonal term\n",level,j,neq);
	    assert( E->BI[level][j] != 0 /* diagonal of matrix = 0, not acceptable */);
            E->BI[level][j]  = (float) 1.0/E->BI[level][j];   
	    }

      if (E->control.verbose)   {
        fprintf(stderr,"output stiffness matrix!!!\n");
        fprintf(E->fp,"level %d\n",level);
        for(j=1;j<=nno;j++)     {
            if(E->NODE[level][j] & OFFSIDE)
                continue;
            eqn1=E->ID[level][j].doff[1];
            eqn2=E->ID[level][j].doff[2];
            max_eqn = E->Node_eqn[level][j+1]-E->Node_eqn[level][j];
            B1=E->Eqn_k[level]+E->Node_k_id[level][j];
            B2=E->Eqn_k[level]+E->Node_k_id[level][j]+max_eqn;
            for(i=0;i<max_eqn;i++) {
                fprintf(E->fp,"%d %d %g %g\n",j,i,B1[i],B2[i]);
                }
            }
        }


    }



    return;
}



/* ============================================
   Function to set up the boundary condition
   masks and other indicators.
   ============================================  */

void construct_masks(E)		/* Add lid/edge masks/nodal weightings */
     struct All_variables *E;
{	
  int i,j,k,l,node,el,elt;
  int lev,elx,elz,ely,nno,nox,noz,noy;
  
  for(lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--){
      elz = E->mesh.ELZ[lev];
      ely = E->mesh.ELY[lev];
      elx = E->mesh.ELX[lev];
      noy = E->mesh.NOY[lev];
      noz = E->mesh.NOZ[lev];
      nno = E->mesh.NNO[lev];

      for(i=1;i<=E->mesh.NNO[lev];i++)
          E->TW[lev][i] = 0.0;

      for(i=1;i<=elz;i++)
          for(j=1;j<=elx;j++)
              for(k=1;k<=ely;k++)
                  for(l=1;l<=enodes[E->mesh.nsd];l++) {
                      elt =  i + (j-1) * elz + (k-1) * elz * elx;
                      node = E->IEN[lev][elt].node[l];
                      E->TW[lev][node] += 1.0;
                  }

      for(i=1;i<=E->mesh.NNO[lev];i++) {
          if(E->NODE[lev][i] & OFFSIDE)
              continue;
          assert( E->TW[lev][i] != 0.0  /* setting weightings failed */);
          E->TW[lev][i] = 1.0/(E->TW[lev][i]);
      }
  }


                                /* Edge masks  */
 

  for(i=1;i<=E->mesh.nox;i++)   /* Horizontal  */
    { for(j=1;j<=E->mesh.noy;j++)
        { node = 1+(i-1)*E->mesh.noz+(j-1)*E->mesh.noz*E->mesh.nox;
          E->node[node] = E->node[node] | TZEDGE;
          E->node[node] = E->node[node] | VZEDGE;
          node += E->mesh.noz-1;;
          E->node[node] = E->node[node] | VZEDGE;
          E->node[node] = E->node[node] | TZEDGE;
        }
    }
  if (E->mesh.nsd == 3) /* not appropriate otherwise */
    for(i=1;i<=E->mesh.noz;i++) /* vertical edge, x normal */
      { for(j=1;j<=E->mesh.noy;j++)
          { node = i + (j-1) * E->mesh.nox * E->mesh.noz;
            E->node[node] = E->node[node] | TXEDGE;
            E->node[node] = E->node[node] | VXEDGE;
            node = i+(E->mesh.nox-1)*E->mesh.noz + (j-1) * E->mesh.nox * E->mesh
.noz;
            E->node[node] = E->node[node] | TXEDGE;
            E->node[node] = E->node[node] | VXEDGE; } }

  for(i=1;i<=E->mesh.noz;i++)   /* vertical edge, y normal */
    { for(j=1;j<=E->mesh.nox;j++)
        { node = i + (j-1) * E->mesh.noz;
          E->node[node] = E->node[node] | TYEDGE;
          E->node[node] = E->node[node] | VYEDGE;
          node = i+(E->mesh.noy-1)*E->mesh.noz*E->mesh.nox + (j-1) * E->mesh.noz
;
          E->node[node] = E->node[node] | TYEDGE;
          E->node[node] = E->node[node] | VYEDGE; } }
 

  return;
  }


/*   ==========================================
     build the sub-element reference matrices
     ==========================================   */

void construct_sub_element(E)
     struct All_variables *E;

{    int i,j,k,l;
     int lev,elx,elz,ely,elzu,elxu,elt,eltu;

     
     for(lev=E->mesh.levmax-1;lev>=E->mesh.levmin;lev--)  {
      elx = E->mesh.ELX[lev];
	  elz = E->mesh.ELZ[lev];
	  ely = E->mesh.ELY[lev];
	  elzu = 2 * elz;
	  elxu = 2 * elx;

	  for(i=1;i<=elx;i++)
	    for(j=1;j<=elz;j++)
	      for(k=1;k<=ely;k++)    {
		  elt = j + (i-1)*elz +(k-1)*elz*elx;
		  eltu = (j*2-1) + elzu *2*(i-1) + elxu*elzu*2*(k-1);

		  for(l=1;l<=enodes[E->mesh.nsd];l++)   {
		      E->EL[lev][elt].sub[l] = eltu
			+ offset[l].vector[1] 
			  + offset[l].vector[0] * elzu
			    + offset[l].vector[2] * elzu * elxu; 
		      }
		  }  
	  }
       
     return; }


void construct_elt_ks(E)
     struct All_variables *E;
{ 
    int e,el,lev,j,k,ii;
    void get_elt_k();
    void matrix_transform_K();
    void get_aug_k();
    void build_diagonal_of_K();

    const int dims=E->mesh.nsd;
    const int n=loc_mat_size[E->mesh.nsd];

    if(E->control.verbose )
	fprintf(stderr,"storing elt k matrices\n");

    for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)  {

	for(el=1;el<=E->mesh.NEL[lev];el++)    {

	    get_elt_k(E,el,E->elt_k[lev][el].k,lev,0);  /* not for penalty */ 

	    if (E->control.augmented_Lagr)
	        get_aug_k(E,el,E->elt_k[lev][el].k,lev,0);

        build_diagonal_of_K(E,el,E->elt_k[lev][el].k,lev);


	    }

            for(j=0;j<E->mesh.NEQ[lev];j++) {
	       if(E->BI[lev][j] ==0.0)  fprintf(stderr,"level %d, equation %d/%d has zero diagonal term\n",lev,j,E->mesh.NEQ[lev]);
               assert( E->BI[lev][j] != 0 /* diagonal of matrix = 0, not acceptable */);
               E->BI[lev][j]  = (float) 1.0/E->BI[lev][j];    
	       }
	}

 if (E->control.verbose) 
    for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)
	for(el=1;el<=E->mesh.NEL[lev];el++)
             for(j=1;j<=enodes[E->mesh.nsd];j++)
		 for(k=1;k<=enodes[E->mesh.nsd];k++) {
		    ii = (j*n+k)*dims-(dims*n+dims);
	/*  fprintf(E->fp,"stiff_for_e %d %d %d %g %g %g %g \n",el,j,k,E->elt_k[lev][el].k[ii],E->elt_k[lev][el].k[ii+1],E->elt_k[lev][el].k[ii+n],E->elt_k[lev][el].k[ii+n+1]);      */
		}

  return; }



void construct_elt_gs(E)
     struct All_variables *E;
{ int el,lev,a;
  void get_elt_g();
  void matrix_transform_g();

  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];

  if(E->control.verbose)
      fprintf(stderr,"storing elt g matrices\n");

  for(lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++)
    for(el=1;el<=E->mesh.NEL[lev];el++)       {
      get_elt_g(E,el,E->elt_del[lev][el].g,lev);  
      }

  return; }



void construct_mat_group(E)
     struct All_variables *E;
{ 

  int llayer,layers(),i,j,k,el,lev,a,nodea,crit2;
  float slope1,slope2,x2;

  const int  dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];
  
  crit2= ends/2;
  for (el=1; el<=E->mesh.nel; el++)  {
    E->mat[el] = 1;
    x2 = 0.0;
    for (a=1;a<=ends;a++)   {
       nodea = E->ien[el].node[a];
       x2 += E->X[2][nodea];
       }
    x2 = x2/ends; 
    llayer = layers(E,x2);
    if (llayer)  {
         E->mat[el]=llayer;
         }
    }

// use cmb values as reference
  

  E->data.therm_exp_factor = 1.0/E->data.therm_exp_factor;

  
  slope1 = (1.0-E->data.therm_exp_factor)/(E->sphere.ro-E->sphere.ri);
  slope2 = (1.0-E->data.therm_diff_factor)/(E->sphere.ro-E->sphere.ri);
  for (i=1;i<=E->mesh.noz;i++)  {
    E->expansivity[i] = (slope1*(E->X[2][i]-E->sphere.ri) + E->data.therm_exp_factor); 
    E->diffusivity[i] = (slope2*(E->X[2][i]-E->sphere.ri) + E->data.therm_diff_factor); 
    }

 for (i=1;i<=E->mesh.noz;i++)
   fprintf(E->fp,"%d  %g %g\n",i,E->expansivity[i],E->diffusivity[i]);


/*
  for (el=1; el<=E->mesh.nel; el++)
    fprintf(E->fp,"mat[%d]= %d \n",el,E->mat[el]);
*/
  return; 
  }


/* routine for constructing stiffness and node_maps */

void construct_stiffness_B_matrix(E)
  struct All_variables *E;
{
  void vcopy();
  void build_diagonal_of_K();
  void build_diagonal_of_Ahat();
  void project_viscosity();
  void construct_node_maps();
  void construct_node_ks(); 
  void construct_elt_ks(); 
  void construct_elt_gs(); 
  
  static int been_here = 0;
  static int been_here0 = 0;

  int i;

  if (been_here0 == 0)  {
   for (i=E->mesh.levmin;i<=E->mesh.levmax;i++)
    if(!E->control.NMULTIGRID && !E->control.NASSEMBLE)  {
       E->elt_k[i]=(struct EK *)malloc((E->mesh.NEL[i]+1)*sizeof(struct EK));
       }
    else if(E->control.NMULTIGRID || E->control.NASSEMBLE) {
       E->Node_eqn[i]  = (int *) malloc((E->mesh.NNO[i]+5) * sizeof(int));
       E->Node_k_id[i] = (int *) malloc((E->mesh.NNO[i]+5) * sizeof(int));
       }
    }
    
  if (been_here0 == 0 || E->viscosity.update_allowed)   {
     /* do the following for the 1st time or update_allowed is true */

    if (E->control.NMULTIGRID)
       project_viscosity(E);
    
    construct_elt_gs(E);
   
    if (E->control.NMULTIGRID || E->control.NASSEMBLE) {
      if (been_here == 0)   {    /* node_maps only built once */
        construct_node_maps(E);
        been_here = 1;
        }
      construct_node_ks(E);
      } 
    else {
      construct_elt_ks(E);
      }
     
    for(i=E->mesh.levmax;i>=E->mesh.levmin;i--)   
      build_diagonal_of_Ahat(E,i);
 
    }

    been_here0 = 1;

  return;
}
