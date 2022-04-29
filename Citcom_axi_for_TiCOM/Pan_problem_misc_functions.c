				
#ifndef __ELEMENT_DEFINITIONS_H__
#define __ELEMENT_DEFINITIONS_H__
#include "element_definitions.h"
#endif

#ifndef __GLOBAL_DEFS_H__
#define __GLOBAL_DEFS_H__
#include "global_defs.h"
#endif

#include <stdio.h>


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
#include <unistd.h>

#if (defined __sunos__)
#define STRING_H
#include <string.h>
#else
#include <strings.h>
#endif

#ifndef STRING_H
#define STRING_H
#include <string.h>
#endif

#if defined(__sgi) || defined(__osf__)
#include <sys/types.h>
#endif


int get_process_identifier()   
{
    int pid;

    pid = (int) getpid();
    return(pid);
}

void unique_copy_file(E,name,comment)
    struct All_variables *E;
    char *name, *comment;
{
    char unique_name[500];
    char command[600];

    sprintf(unique_name,"%s-%06d.%s",name,E->control.PID,comment);
    sprintf(command,"cp -f \"%s\" \"%s\"\n",name,unique_name);
    system(command);
}


void thermal_buoyancy(E)
     struct All_variables *E;
    
{ 
    int i,j;
    float *H,slope,temp1;
    float g,r,Rc,Ddens;
    float buoyancy_enhancement; /* due to different material props in l.m. */
    void remove_horiz_ave();
    void phase_change();

    H = (float *)malloc( (E->mesh.noz+1)*sizeof(float));

    buoyancy_enhancement = E->data.density_lm/E->data.density;

    if (abs(E->control.Ra_670)>0.0 || abs(E->control.Ra_410)>0.0)   {

      phase_change(E,E->Fas670,E->Fas670_b,E->Fas410,E->Fas410_b);

      for(i=1;i<=E->mesh.nno;i++)
        E->buoyancy[i] = -E->control.Ra_670 * E->Fas670[i] - E->control.Ra_410 * E->Fas410[i];
      }
    else
      for(i=1;i<=E->mesh.nno;i++)
        E->buoyancy[i] = 0.0;

   for(i=1;i<=E->mesh.nno;i++)   {

/*  for constant density planet with g decreasing linearly with depth */
/*        E->buoyancy[i] = E->control.Ra_temp * E->T[i] * E->X[2][i]
                       + E->control.Ra_comp * E->C[i] * E->X[2][i];
*/

/*  for constant g in the mantle */

       j = (i-1)%(E->mesh.noz)+1;
       E->buoyancy[i] += E->control.Ra_temp * E->T[i] 
                      * E->expansivity[j]
                      - E->control.Ra_comp * E->C[i]; 

/*  buoyancy in silicate portion different from that in ice due to 
 *  different material properties.  Change that here.
 *  This may not be the best way to do this. Better to use E->C ?
 */

	if(E->X[2][i] <= E->viscosity.nlm)
		E->buoyancy[i] *= buoyancy_enhancement;

/* for g(r) */
/*
       j = (i-1)%(E->mesh.noz)+1;
	Ddens = E->data.density_core/E->data.density - 1.0;
	Rc = E->sphere.ri;
	r = E->X[2][i];
	g = ((Ddens*Rc*Rc*Rc)/(r*r) + r) / ((Ddens*Rc*Rc*Rc) + 1.0);
	E->buoyancy[i] += (E->control.Ra_temp * E->T[i] * E->expansivity[j]
                      - E->control.Ra_comp * E->C[i]) * g;

	if(i<=E->mesh.noz && E->monitor.solution_cycles==0)
		fprintf(stderr,"%d %f %f %f %f\n",i,r,g,Ddens,Rc);
*/

  /* Additional buoyancy due to isostatic uplift in transient crater */
  /* Apply ONLY when an impact occurs during the timestep */
  if (E->impacts.now >= 0) {
    E->buoyancy[i] += 1.0*E->control.Ra_temp * (E->impacts.H_t[i]/E->sphere.ro_dim) / (E->sphere.ro - E->sphere.ri);
    if(i==1) {
      fprintf(stderr,"Isostatic rebound!!!\n");
      fprintf(E->fp,"Isostatic rebound!!!\n");
      }
    }

  }

    remove_horiz_ave(E,E->buoyancy,H,0);

  free ((void *) H);
    return;
}
 
double SIN_D(x)
     double x;
{
#if defined(__osf__)
  return sind(x);
#else
  return sin((x/180.0) * M_PI);
#endif

}

double COT_D(x)
     double x;
{
#if defined(__osf__)
  return cotd(x);
#else
  return tan(((90.0-x)/180.0) * M_PI);
#endif

}


/* non-runaway malloc */

void * Malloc1(bytes,file,line)
    int bytes;
    char *file;
    int line;
{
    void *ptr;

    ptr = malloc((size_t)bytes);
    if (ptr == (void *)NULL) {
	fprintf(stderr,"Memory: cannot allocate another %d bytes \n(line %d of file %s)\n",bytes,line,file);
	exit(0);
    }

    return(ptr);
}


/* Read in a file containing previous values of a field. The input in the parameter
   file for this should look like: `previous_name_file=string' and `previous_name_column=int' 
   where `name' is substituted by the argument of the function. 

   The file should have the standard CITCOM output format:
     # HEADER LINES etc
     index X Z Y ... field_value1 ...
     index X Z Y ... field_value2 ...
   where index is the node number, X Z Y are the coordinates and
   the field value is in the column specified by the abbr term in the function argument

   If the number of nodes OR the XZY coordinates for the node number (to within a small tolerance)
   are not in agreement with the existing mesh, the data is interpolated. 

   */

int read_previous_field(E,field,name,abbr)
    struct All_variables *E;
    float *field;
    char *name, *abbr;
{
    void fcopy_interpolating();
    int input_string();

    float cross2d();
    void vcopy();

    char discard[5001];
    char *token;
    char *filename;
    char *input_token;
    FILE *fp;
    int fnodesx,fnodesz,fnodesy;
    int i,j,column,found;
    int interpolate=0;

    float *X,*Z,*Y,*T;

    filename=(char *)malloc(500*sizeof(char));
    input_token=(char *)malloc(1000*sizeof(char));

    /* Define field name, read parameter file to determine file name and column number */

    sprintf(input_token,"previous_%s_file",name);
    if(!input_string(input_token,filename,"initialize")) {
	fprintf(E->fp,"No previous %s information found in input file\n",name);fflush(E->fp);
	return(0);   /* if not found, take no further action, return zero */
    }

     
    fprintf(E->fp,"Previous %s information is in file %s\n",name,filename);fflush(E->fp);
 
    /* Try opening the file, fatal if this fails too */

    if((fp=fopen(filename,"r")) == NULL) {
	fprintf(E->fp,"Unable to open the required file `%s' (this is fatal)",filename);fflush(E->fp);
	if(E->control.verbose)
	   	fprintf(stderr,"Unable to open the required file `%s'",filename);
	exit(1);
    }
  
    
     /* Read header, get nodes xzy */

    fgets(discard,4999,fp);
    fgets(discard,4999,fp); 
    i=sscanf(discard,"# NODESX=%d NODESZ=%d NODESY=%d",&fnodesx,&fnodesz,&fnodesy);
    if(i<3) {
	fprintf(E->fp,"File %s is not in the correct format\n",filename);fflush(E->fp);
	exit(1);
    }

    fgets(discard,4999,fp); /* largely irrelevant line */
    fgets(discard,4999,fp);
    
    /* this last line is the column headers, we need to search for the occurence of abbr to
       find out the column to be read in */

    if(strtok(discard,"|")==NULL) { 
	fprintf(E->fp,"Unable to deciphre the columns in the input file");fflush(E->fp);
	exit(1);
    }

    found=0;
    column=1;

    while(found==0 && (token=strtok(NULL,"|")) != NULL) {
	if(strstr(token,abbr)!=0)
	    found=1;
	column++;
    }

    if(found) {
	fprintf(E->fp,"\t%s (%s) found in column %d\n",name,abbr,column);fflush(E->fp);
    }    
    else {
	fprintf(E->fp,"\t%s (%s) not found in file: %s\n",name,abbr,filename);fflush(E->fp);
	exit(1);
    }
    

  
    /* Another fatal condition (not suitable for interpolation: */
    if(((3!= E->mesh.nsd) && (fnodesy !=1)) || ((3==E->mesh.nsd) && (1==fnodesy))) {
	fprintf(E->fp,"Input data for file `%s'  is of inappropriate dimension (not %dD)\n",filename,E->mesh.nsd);fflush(E->fp);
	exit(1);
    }

    X=(float *)malloc((2+fnodesx*fnodesz*fnodesy)*sizeof(float));
    Z=(float *)malloc((2+fnodesx*fnodesz*fnodesy)*sizeof(float));
    Y=(float *)malloc((2+fnodesx*fnodesz*fnodesy)*sizeof(float));
    T=(float *)malloc((2+fnodesx*fnodesz*fnodesy)*sizeof(float));
    
   /* Format for reading the input file (including coordinates) */

    sprintf(input_token," %%d %%e %%e %%e");
    for(i=5;i<column;i++)
	strcat(input_token," %*f");
    strcat(input_token," %f");


    for(i=1;i<=fnodesx*fnodesz*fnodesy;i++) {
	fgets(discard,4999,fp);
	sscanf(discard,input_token,&j,&(X[i]),&(Z[i]),&(Y[i]),&T[i]);
    }
    /* check consistency & need for interpolation */

    fclose(fp);

    if(fnodesx != E->mesh.nox || fnodesz != E->mesh.noz || fnodesy != E->mesh.noy)
	interpolate=1;
	
    for(i=1;i<=fnodesx*fnodesz*fnodesy;i++)
	if( fabs(X[i]-E->X[1][i]) > 0.01*fabs(X[i]) ||
	    fabs(Z[i]-E->X[2][i]) > 0.01*fabs(Z[i]) ||
	    ((3==E->mesh.nsd) && fabs(Y[i]-E->X[3][i]) > 0.01*fabs(Y[i]))) {
	    interpolate=i;
	    break;
	}

    if(interpolate!=0) {
	fprintf(E->fp,"\t%s requires interpolation from previous value\n",name,interpolate);fflush(E->fp);
	fprintf(E->fp,"\tOld nodes = %d/%d/%d and new nodes = %d/%d/%d\n",fnodesx,fnodesz,fnodesy,E->mesh.nox,E->mesh.noz,E->mesh.noy);fflush(E->fp);
	fcopy_interpolating(E,X,Z,Y,fnodesx,fnodesz,fnodesy,T,field);
    }
    else {
	fprintf(E->fp,"\t%s requires no interpolation from previous value\n",name);fflush(E->fp);
	vcopy(field,T,1,E->mesh.nno);
    }

    free((void *)X);
    free((void *)Z);
    free((void *)Y);
    free((void *)T);
    free((void *)filename);
    free((void *)input_token);
    
    return(1);
}

/* Copy one field to another on a different (but similarly structured) mesh.  The
   field TT (output) is on the standard mesh for this problem, T (input) is on the different
   mesh whose coordinates are described by X,Z,Y 


   */



void fcopy_interpolating(E,X,Z,Y,nx,nz,ny,T,TT)
    struct All_variables *E;
    float *X,*Z,*Y;
    int nx,nz,ny;
    float *T,*TT;
{
    float cross2d();
    void p_to_nodes();
    void p_to_centres();
    void vcopy();
    double CPU_time0(),time;
    
    float *P;

    int i,j,found,not_found;
    int elX,elY,elZ,ex,ez;
    int old_ex,old_ez,old_ey;
    int node1,node2,node3,node4;
    int node5,node6,node7,node8;
    float inside1,inside2,inside3,inside4;
    float inside5,inside6,inside7,inside8,inside9,inside10,inside11,inside12;
    float distance1,distance2,distance3,distance4;
    float distance5,distance6,distance7,distance8;
    float d1,d2,d3,d4,d5,d6,d7,d8;

    const int dims=E->mesh.nsd;
 
    P=(float *)malloc((E->mesh.nno+1)*sizeof(float));
        
    /* run over all the data points (take care to hit only one element), determine
       inside/outside of each element. The criterion for inside-ness is that the
       cross product of the vectors joining the point to the ends of each edge should
       have the same sign as the dot product of the centre to the ends of each edge.
       There are, undoubtedly, better ways to do this !!!

       Because a CITCOM node-ordering is assumed, we can also guess that the best place to start
       looking for a node is where you found the last one !

       */

    old_ex=old_ez=old_ey=1;

    elX = nx-1;
    elZ = nz-1;


     if(E->control.print_convergence) {   
	 time=CPU_time0();
	 fprintf(stderr,"Interpolating ...");
     }

    not_found=0;
  
    if(2==dims)
	for(i=1;i<=E->mesh.nno;i++) {
	    found=0;
	    for(ex=old_ex;ex<=elX && found==0;ex++)
		for(ez=1;ez<=elZ && found==0 ;ez++) {
		    node1=ez+(ex-1)*nz;
		    node2=node1+offset[2].vector[1]+offset[2].vector[0]*nz;
		    node3=node1+offset[3].vector[1]+offset[3].vector[0]*nz;
		    node4=node1+offset[4].vector[1]+offset[4].vector[0]*nz;
		    
		    if ((inside1 = cross2d(X[node1]-E->X[1][i],Z[node1]-E->X[2][i],X[node2]-E->X[1][i],Z[node2]-E->X[2][i],3)) <= 0.0 &&
			(inside4 = cross2d(X[node4]-E->X[1][i],Z[node4]-E->X[2][i],X[node1]-E->X[1][i],Z[node1]-E->X[2][i],3)) <= 0.0 &&
			(inside2 = cross2d(X[node2]-E->X[1][i],Z[node2]-E->X[2][i],X[node3]-E->X[1][i],Z[node3]-E->X[2][i],3)) <= 0.0 &&
			(inside3 = cross2d(X[node3]-E->X[1][i],Z[node3]-E->X[2][i],X[node4]-E->X[1][i],Z[node4]-E->X[2][i],3)) <= 0.0) {
			found = node1;
			old_ex=ex;
		    }
			
		}
	    
	    /* finish the loop if not found */
	    for(ex=1;ex<=old_ex && found==0;ex++)
		for(ez=1;ez<=elZ && found==0 ;ez++) {
		    node1=ez+(ex-1)*nz;
		    node2=node1+offset[2].vector[1]+offset[2].vector[0]*nz;
		    node3=node1+offset[3].vector[1]+offset[3].vector[0]*nz;
		    node4=node1+offset[4].vector[1]+offset[4].vector[0]*nz;
		    
		    if ((inside1 = cross2d(X[node1]-E->X[1][i],Z[node1]-E->X[2][i],X[node2]-E->X[1][i],Z[node2]-E->X[2][i],3)) <= 0.0 &&
			(inside4 = cross2d(X[node4]-E->X[1][i],Z[node4]-E->X[2][i],X[node1]-E->X[1][i],Z[node1]-E->X[2][i],3)) <= 0.0 &&
			(inside2 = cross2d(X[node2]-E->X[1][i],Z[node2]-E->X[2][i],X[node3]-E->X[1][i],Z[node3]-E->X[2][i],3)) <= 0.0 &&
			(inside3 = cross2d(X[node3]-E->X[1][i],Z[node3]-E->X[2][i],X[node4]-E->X[1][i],Z[node4]-E->X[2][i],3)) <= 0.0) {
			found = node1;
			old_ex=ex;
		    }
		}

	    /* and having found the right node location, interpolate the appropriate value to it */

	    if(!found)
		not_found++;
	    else {
		distance1 = ((X[node1]-E->X[1][i])*(X[node1]-E->X[1][i])+(Z[node1]-E->X[2][i])*(Z[node1]-E->X[2][i]));
		distance2 = ((X[node2]-E->X[1][i])*(X[node2]-E->X[1][i])+(Z[node2]-E->X[2][i])*(Z[node2]-E->X[2][i]));
		distance3 = ((X[node3]-E->X[1][i])*(X[node3]-E->X[1][i])+(Z[node3]-E->X[2][i])*(Z[node3]-E->X[2][i]));
		distance4 = ((X[node4]-E->X[1][i])*(X[node4]-E->X[1][i])+(Z[node4]-E->X[2][i])*(Z[node4]-E->X[2][i]));
	    
		d1=distance2*distance3*distance4;
		d2=distance1*distance3*distance4;
		d3=distance2*distance1*distance4;
		d4=distance2*distance3*distance1;
		
		TT[i] = (d1*T[node1]+d2*T[node2]+d3*T[node3]+d4*T[node4])/(d1+d2+d3+d4);
	    
	    }
	}

    else {
	elY = (3==dims)? ny-1 : 1;
	fprintf(stderr,"3D interolator not yet implemented !!! \n");
   	vcopy(TT,T,1,E->mesh.nno); }
   
    if(E->control.print_convergence)
	fprintf(stderr,". done (%f secs)\n",CPU_time0()-time);

    if(not_found)
	fprintf(E->fp,"Warning: unable to interpolate old  data to %d nodes in the new mesh\n",not_found);
    else {
	p_to_centres(E,TT,P,E->mesh.levmax); /* if interpolated, apply slight smoothing */
	p_to_nodes(E,P,TT,E->mesh.levmax);
    } 

    free((void *)P);

    return;
}


/* returns the out of plane component of the cross product of
   the two vectors assuming that one is looking AGAINST the 
   direction of the axis of D, anti-clockwise angles
   are positive (are you sure ?), and the axes are ordered 2,3 or 1,3 or 1,2 */


float cross2d(x11,x12,x21,x22,D)
    float x11,x12,x21,x22;
    int D;
{
   if(1==D)
       return( x11*x22-x12*x21);
   if(2==D) 
       return(-x11*x22+x12*x21);
   if(3==D)
       return( x11*x22-x12*x21);

}


void field_arbitrary_rectangle_file(E,parse_and_apply,RECT,name,field,BC,bcbitf,bcmask_on,bcmask_off)
    struct All_variables *E;
    int parse_and_apply;
    struct Rect *RECT;
    char *name;
    float *field;
    int BC;
    unsigned int *bcbitf;
    unsigned int bcmask_on,bcmask_off;
{
    
    char read_string[500];
    float weight,radius2,weight2;
    float x1,y1,z1;
    int in1,in2,in3;
    int number,node;
    int combine_option;

    void field_arbitrary_rectangle();
    int input_int();
    int input_float_vector();
    int input_char_vector();

    sprintf(read_string,"%s_rect",name);
    input_int(read_string,&(RECT->numb),"0");
    sprintf(read_string,"%s_rectx1",name);
    input_float_vector(read_string,RECT->numb,RECT->x1);
    sprintf(read_string,"%s_rectx2",name);
    input_float_vector(read_string,RECT->numb,RECT->x2);
    sprintf(read_string,"%s_rectz1",name);
    input_float_vector(read_string,RECT->numb,RECT->z1);
    sprintf(read_string,"%s_rectz2",name);
    input_float_vector(read_string,RECT->numb,RECT->z2);
    sprintf(read_string,"%s_recty1",name);
    input_float_vector(read_string,RECT->numb,RECT->y1);
    sprintf(read_string,"%s_recty2",name);
    input_float_vector(read_string,RECT->numb,RECT->y2);
    sprintf(read_string,"%s_recthw",name);
    input_float_vector(read_string,RECT->numb,RECT->halfw);
    sprintf(read_string,"%s_rectmag",name);
    input_float_vector(read_string,RECT->numb,RECT->mag);
    sprintf(read_string,"%s_rectovl",name);
    input_char_vector(read_string,RECT->numb,RECT->overlay);
 
    if(parse_and_apply)
	field_arbitrary_rectangle(E,RECT,field,BC,bcbitf,bcmask_on,bcmask_off);
 
    return;  
}

void field_arbitrary_rectangle(E,RECT,field,BC,bcbitf,bcmask_on,bcmask_off)
    struct All_variables *E;
    struct Rect *RECT;
    float *field;
    int BC;
    unsigned int *bcbitf;
    unsigned int bcmask_on,bcmask_off;
{
    
    float weight,radius2,weight2;
    float x1,y1,z1;
    int in1,in2,in3;
    int number,node;
    int combine_option;
  
    for(node=1;node<=E->mesh.nno;node++) {
	x1=E->X[1][node];
	z1=E->X[2][node];
	y1=(E->mesh.nsd!=3) ? 0.0 : E->X[3][node];
	
	for(number=0;number<RECT->numb;number++) {
	   
	    switch (RECT->overlay[number]) {
	    case 'M':
		weight=1.0;
		break;
	    case 'R':
	    case 'A':
		weight=0.0;
		break;
	    }

	    in1=(x1 >= RECT->x1[number] && x1 <= RECT->x2[number]);
	    in2=(z1 >= RECT->z1[number] && z1 <= RECT->z2[number]);
	    in3=(3!=E->mesh.nsd || y1 >= RECT->y1[number] && y1 <= RECT->y2[number]);
	    
	    if(in1 && in2 && in3) {
		weight = RECT->mag[number];
		radius2=0.0;
	    }
	    else {
		radius2 = 
		    (in1 ? 0.0 : min(1.0e-10+(x1-RECT->x1[number])*(x1-RECT->x1[number]),
				     1.0e-10+(x1-RECT->x2[number])*(x1-RECT->x2[number]))) +
		    (in2 ? 0.0 : min(1.0e-10+(z1-RECT->z1[number])*(z1-RECT->z1[number]),
				     1.0e-10+(z1-RECT->z2[number])*(z1-RECT->z2[number]))) +
		    (in3 ? 0.0 : min(1.0e-10+(y1-RECT->y1[number])*(y1-RECT->y1[number]),
				     1.0e-10+(y1-RECT->y2[number])*(y1-RECT->y2[number])));
		        
		weight += RECT->mag[number]*exp(-radius2/(0.5 * RECT->halfw[number]*RECT->halfw[number]));
	    }
		
	    switch (RECT->overlay[number]) {
	    case 'R':
		if(radius2 > RECT->halfw[number]* RECT->halfw[number])
		    break;
		weight2=1.0/(1.0+radius2/(1.0e-10+RECT->halfw[number]* RECT->halfw[number]));
		field[node]=(1.0-weight2)*field[node]+weight2*weight;
		/*fprintf(stderr," %d (%f/%f) = %f,%f   /   %f (%f)\n",node,x1,z1,weight,weight2,field[node],radius2);*/
		break;
	    case 'M':
		field[node]*=weight;
		break;
	    case 'A':
		field[node]+=weight;
		break;
	    default:
		fprintf(E->fp,"RECTANGLE: %d can't work out how to combine new/old fields\n",number);
		break;
	    }
		
	    
	    if(BC) {
		bcbitf[node] = (bcbitf[node] | bcmask_on);
		bcbitf[node] = (bcbitf[node] & (~bcmask_off));
	    }
	}

   }


    return;  
}


void field_arbitrary_circle_file(E,parse_and_apply,CIRC,name,field,BC,bcbitf,bcmask_on,bcmask_off)
    struct All_variables *E;
    int parse_and_apply;
    struct Circ *CIRC;
    char *name;
    float *field;
    int BC;
    unsigned int *bcbitf;
    unsigned int bcmask_on,bcmask_off;
{
    
    char read_string[500];
    float weight,radius2,weight2;
    float x1,y1,z1;
    int in1;
    int number,node;
    int combine_option;
    void field_arbitrary_circle();
    int input_int();
    int input_float_vector();
    int input_char_vector();
    
    sprintf(read_string,"%s_circ",name);
    input_int(read_string,&(CIRC->numb),"0");
    sprintf(read_string,"%s_circx",name);
    input_float_vector(read_string,CIRC->numb,CIRC->x);
    sprintf(read_string,"%s_circz",name);
    input_float_vector(read_string,CIRC->numb,CIRC->z);
    sprintf(read_string,"%s_circy",name);
    input_float_vector(read_string,CIRC->numb,CIRC->y);
    sprintf(read_string,"%s_circrad",name);
    input_float_vector(read_string,CIRC->numb,CIRC->rad);
    sprintf(read_string,"%s_circmag",name);
    input_float_vector(read_string,CIRC->numb,CIRC->mag);
    sprintf(read_string,"%s_circhw",name);
    input_float_vector(read_string,CIRC->numb,CIRC->halfw);
    sprintf(read_string,"%s_circovl",name);
    input_char_vector(read_string,CIRC->numb,CIRC->overlay);
 
    if(parse_and_apply)
	field_arbitrary_circle(E,CIRC,field,BC,bcbitf,bcmask_on,bcmask_off);
   
    return;  
}
void field_arbitrary_circle(E,CIRC,field,BC,bcbitf,bcmask_on,bcmask_off)
    struct All_variables *E;
    struct Circ *CIRC;
    float *field;
    int BC;
    unsigned int *bcbitf;
    unsigned int bcmask_on,bcmask_off;
{
    
    char read_string[500];
    float weight,radius2,weight2;
    float x1,y1,z1;
    int in1;
    int number,node;
    int combine_option;
  
   
    for(node=1;node<=E->mesh.nno;node++) {
	x1=E->X[1][node];
	z1=E->X[2][node];
	y1=(E->mesh.nsd!=3) ? 0.0 : E->X[3][node];
	
	for(number=0;number<CIRC->numb;number++) {
	    switch (CIRC->overlay[number]) {
	    case 'M':
		weight=1.0;
		break;
	    case 'R':
	    case 'A':
		weight=0.0;
		break;
	    }
	    
	    radius2 =
		(x1-CIRC->x[number])*(x1-CIRC->x[number]) +
		(z1-CIRC->z[number])*(z1-CIRC->z[number]) +
		((E->mesh.nsd != 3) ? 0.0 : (y1-CIRC->y[number])*(y1-CIRC->y[number]));
 
	     if(radius2 <= CIRC->rad[number]*CIRC->rad[number]) {
		weight = CIRC->mag[number];
		radius2=0.0;
	    }
	    else {
		radius2 -= CIRC->rad[number] * CIRC->rad[number];
		weight += CIRC->mag[number]*exp(-2.0*radius2/(1.0e-10+CIRC->halfw[number]*CIRC->halfw[number]));
	    }
		
	    switch (CIRC->overlay[number]) {
	    case 'R':
		if(radius2 > CIRC->halfw[number]*CIRC->halfw[number])
		    break;
		weight2=1.0/(1.0+radius2/(1.0e-10+CIRC->halfw[number]*CIRC->halfw[number]));
		field[node]=(1.0-weight2)*field[node]+weight2*weight;

		/*fprintf(stderr," %d (%f/%f) = %f,%f   /   %f (%f)\n",node,x1,z1,weight,weight2,field[node],radius2);*/
		break;
	    case 'M':
		field[node]*=weight;
		break;
	    case 'A':
		field[node]+=weight;
		break;
	    default:
		fprintf(E->fp,"CIRCLE: %d can't work out how to combine new/old fields\n",number);
		break;
	    }
		
	    
	    if(BC) {
		bcbitf[node] = (bcbitf[node] | bcmask_on);
		bcbitf[node] = (bcbitf[node] & (~bcmask_off));
	    }
	}

   }


    return;  
}

		
void field_arbitrary_harmonic_file(E,parse_and_apply,HARM,name,field,BC,bcbitf,bcmask_on,bcmask_off)  
    struct All_variables *E;
    int parse_and_apply;
    struct Harm *HARM;
    char *name;
    int BC;
    float *field;
    unsigned int *bcbitf;
    unsigned int bcmask_on,bcmask_off;
{
  char read_string[500];
  void field_arbitrary_harmonic();
  int input_int();
  int input_float_vector();
  int input_char_vector();
  int i;

  sprintf(read_string,"%s_harm",name);
  input_int(read_string,&(HARM->numb),"0");
  sprintf(read_string,"%s_harms",name);
  input_int(read_string,&(HARM->harms),"0,0,19");
  sprintf(read_string,"%s_harmoff",name);
  input_float_vector(read_string,HARM->numb,HARM->off);
  sprintf(read_string,"%s_harmx1",name);
  input_float_vector(read_string,HARM->numb,HARM->x1);
  sprintf(read_string,"%s_harmx2",name);
  input_float_vector(read_string,HARM->numb,HARM->x2);
  sprintf(read_string,"%s_harmz1",name);
  input_float_vector(read_string,HARM->numb,HARM->z1);
  sprintf(read_string,"%s_harmz2",name);
  input_float_vector(read_string,HARM->numb,HARM->z2);
  sprintf(read_string,"%s_harmy1",name);
  input_float_vector(read_string,HARM->numb,HARM->z1);
  sprintf(read_string,"%s_harmy2",name);
  input_float_vector(read_string,HARM->numb,HARM->z2);
  sprintf(read_string,"%s_harmovl",name);
  input_char_vector(read_string,HARM->numb,HARM->overlay);
 
  for(i=0;i<HARM->harms;i++) {
      sprintf(read_string,"%s_harmkx%02d",name,i+1);
      input_float_vector(read_string,HARM->numb,HARM->kx[i]);
      sprintf(read_string,"%s_harmkz%02d",name,i+1);
      input_float_vector(read_string,HARM->numb,HARM->kz[i]);
      sprintf(read_string,"%s_harmky%02d",name,i+1);
      input_float_vector(read_string,HARM->numb,HARM->ky[i]);
      sprintf(read_string,"%s_harmka%02d",name,i+1);
      input_float_vector(read_string,HARM->numb,HARM->ka[i]);
      sprintf(read_string,"%s_harmphx%02d",name,i+1);
      input_float_vector(read_string,HARM->numb,HARM->phx[i]);
      sprintf(read_string,"%s_harmphz%02d",name,i+1);
      input_float_vector(read_string,HARM->numb,HARM->phz[i]);
      sprintf(read_string,"%s_harmphy%02d",name,i+1);
      input_float_vector(read_string,HARM->numb,HARM->phy[i]);
  }
 
  if(parse_and_apply)
	field_arbitrary_harmonic(E,HARM,field,BC,bcbitf,bcmask_on,bcmask_off);
 
  return;  
}

void field_arbitrary_harmonic(E,HARM,field,BC,bcbitf,bcmask_on,bcmask_off)
    struct All_variables *E;
    struct Harm *HARM;
    float *field;
    int BC;
    unsigned int *bcbitf;
    unsigned int bcmask_on,bcmask_off;
{
    
    float weight,radius2,weight2;
    float x1,y1,z1;
    int in1,in2,in3;
    int number,node,l;
    int combine_option;
  
    for(node=1;node<=E->mesh.nno;node++) {
	x1=E->X[1][node];
	z1=E->X[2][node];
	y1=(E->mesh.nsd!=3) ? 0.0 : E->X[3][node];
	
	for(number=0;number<HARM->numb;number++) {
	   
	    switch (HARM->overlay[number]) {
	    case 'M':
		weight=1.0;
		break;
	    case 'R':
	    case 'A':
		weight=0.0;
		break;
	    }

	    in1=(x1 >= HARM->x1[number] && x1 <= HARM->x2[number]);
	    in2=(z1 >= HARM->z1[number] && z1 <= HARM->z2[number]);
	    in3=(3!=E->mesh.nsd || y1 >= HARM->y1[number] && y1 <= HARM->y2[number]);
	    
	    if(in1 && in2 && in3) {
		weight = HARM->off[number];
		for(l=0;l<HARM->harms;l++) {
		    weight += HARM->ka[l][number] *
			cos((HARM->kx[l][number]*x1+HARM->phx[l][number])*M_PI) *
			cos((HARM->kz[l][number]*z1+HARM->phz[l][number])*M_PI) *
			cos((HARM->ky[l][number]*y1+HARM->phy[l][number])*M_PI) ;
		}
	    
	   	
		switch (HARM->overlay[number]) {
		case 'R':
		    field[node]=weight;
		    break;
		case 'M':
		    field[node]*=weight;
		    break;
		case 'A':
		    field[node]+=weight;
		    break;
		default:
		    fprintf(E->fp,"POLYNOMIAL: %d can't work out how to combine new/old fields\n",number);
		    break;
		}
		
	    
		if(BC) {
		    bcbitf[node] = (bcbitf[node] | bcmask_on);
		    bcbitf[node] = (bcbitf[node] & (~bcmask_off));
		}
	    }
	}
    }


    return;  
}


/* =====================================*/
 double sphere_h(l,m,t,f,ic)
 int l,m,ic;
 double t,f;
 {
 
 double plgndr_a(),sphere_hamonics;

 sphere_hamonics = 0.0;
 if (ic==0)
    sphere_hamonics = cos(m*f)*plgndr_a(l,m,t);
 else if (m)
    sphere_hamonics = sin(m*f)*plgndr_a(l,m,t);

 return sphere_hamonics;
 }

/* =====================================*/
 double plgndr_a(l,m,t)
 int l,m;
 double t;
 {

  int i,ll;
  double x,fact,pll,pmm,pmmp1,somx2,plgndr;
  const double two=2.0;
  const double one=1.0;

  x = cos(t);
  pmm=one; 
  if(m>0) {
    somx2=sqrt((one-x)*(one+x));
    fact = one;
    for (i=1;i<=m;i++)   {
      pmm = -pmm*fact*somx2;
      fact = fact + two;
      }
    }

  if (l==m) 
     plgndr = pmm;
  else  {
     pmmp1 = x*(2*m+1)*pmm;
     if(l==m+1)  
       plgndr = pmmp1;
     else   {
       for (ll=m+2;ll<=l;ll++)  {
         pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
         pmm = pmmp1;
         pmmp1 = pll;
         }
       plgndr = pll;
       }
     }

 return plgndr;
 }

/* =====================================
 =====================================*/
 double modified_plgndr_a(l,m,t)
 int l,m;
 double t;
 {

  int i,ll;
  double x,fact1,fact2,fact,pll,pmm,pmmp1,somx2,plgndr;
  const double three=3.0;
  const double two=2.0;
  const double one=1.0;

  x = cos(t);
  pmm=one; 
  if(m>0) {
    somx2=sqrt((one-x)*(one+x));
    fact1= three;
    fact2= two;
    for (i=1;i<=m;i++)   {
      fact=sqrt(fact1/fact2);
      pmm = -pmm*fact*somx2;
      fact1+=  two;
      fact2+=  two;
      }
    }

  if (l==m) 
     plgndr = pmm;
  else  {
     pmmp1 = x*sqrt(two*m+three)*pmm;
     if(l==m+1)  
       plgndr = pmmp1;
     else   {
       for (ll=m+2;ll<=l;ll++)  {
	 fact1= sqrt((4.0*ll*ll-one)*(double)(ll-m)/(double)(ll+m));
	 fact2= sqrt((2.0*ll+one)*(ll-m)*(ll+m-one)*(ll-m-one)
		     /(double)((two*ll-three)*(ll+m)));
         pll = ( x*fact1*pmmp1-fact2*pmm)/(ll-m);
         pmm = pmmp1;
         pmmp1 = pll;
         }
       plgndr = pll;
       }
     }

 plgndr /= sqrt(4.0*M_PI);

 if (m!=0) plgndr *= sqrt(two);

 return plgndr;
 }

 /* ===================================  */
  double sqrt_multis(jj,ii)
  int ii,jj;
 {
  int i;
  double sqrt_multisa;

  sqrt_multisa = 1.0;
  if(jj>ii)
    for (i=jj;i>ii;i--)
      sqrt_multisa *= 1.0/sqrt((double)i);

  return sqrt_multisa;
  }

 /* ===================================  */
  double multis(ii)
  int ii;
 {
  int i;
  double multisa;

  multisa = 1.0;
  if (ii)
    for (i=2;i<=ii;i++)
      multisa *= (double)i;

  return multisa;
  }


 /* ===================================  */
 int int_multis(ii)
 int ii;
 {
 int i,multisa;

 multisa = 1;
 if (ii)
   for (i=2;i<=ii;i++)
     multisa *= i;

 return multisa;
 }
