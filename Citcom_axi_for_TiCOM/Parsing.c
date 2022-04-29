/* Routines which read filenames from the command line and
   then parse the contents as parameters for citcom */

#include "Parsing.h"

#ifndef __EXT_PAR__
#define __EXT_PAR__
struct ext_par		/* global variables for getpar */
{
  char *progname;
  int argflags;
  struct arglist *arglist;
  struct arglist *arghead;
  char *argbuf;
  int nlist;
  int nbuf;
  int listmax;
  int bufmax;
  FILE *listout;
}	ext_par;
#endif

#ifndef __ARGLIST__
#define __ARGLIST__
struct arglist		/* structure of list set up by setpar */
{
    int argname_offset;
    int argval_offset;
    int hash;
};
#endif

#ifndef __VERBOSE_DESCRIBE_BEGINNER__
#define __VERBOSE_DESCRIBE_BEGINNER__
int VERBOSE = 0;
int DESCRIBE = 0;
int BEGINNER = 0;
#endif

#ifndef __SETUP_PARSER__
#define __SETUP_PARSER__
void setup_parser(E,ac,av)
     struct All_variables *E;
     int ac; 
     char **av;     
{
    void unique_copy_file();
    
    FILE * fp;
    char *pl,*pn,*pv;
    char t1, t2, line[MAXLINE], name[MAXNAME], value[MAXVALUE];
    int i,j,k;
  
  
    /* should get file length & cpp &c before any further parsing */

    /* for now, read one filename from the command line, we'll parse that ! */
    
    if (ac < 2)   {
	fprintf(stderr,"Usage: citcom PARAMETERFILE\n");
	exit(10);
    }

  
    if ((fp = fopen(av[1],"r")) == NULL)  {
      fprintf(stderr,"File: %s is unreadable\n",av[1]);
      exit(11);
    }


    unique_copy_file(E,av[1],"copy");

  /* now the parameter file is open, read into memory */
  int add_to_parameter_list();
  while( fgets(line,MAXLINE,fp) != NULL )
    { pl= line;
      /* loop over entries on each line */
    loop:	
      while(*pl==' ' || *pl=='\t') pl++;
      if(*pl=='\0'|| *pl=='\n') continue; /* end of line */
      if(*pl=='#') continue; /* end of interpretable part of line */

      /* get name */
      pn= name;
      while(*pl != '=' && *pl != '\0' && *pl != ' '
	    && *pl != '\n'		/* FIX by Glenn Nelson */
	    && *pl != '\t') 
	*pn++ = *pl++;
      *pn = '\0';
      if(*pl == '=') pl++;
      
      /* get value */
      *value= '\0';
      pv= value;
      if(*pl=='"' || *pl=='\'')
	t1= t2= *pl++; 
      else
	{ t1= ' ';
	  t2= '\t';
	}
      while(*pl!=t1 && *pl!=t2 &&
	    *pl!='\0' && *pl!='\n') *pv++= *pl++;
      *pv= '\0';
      if(*pl=='"' || *pl=='\'')
	pl++;
      add_to_parameter_list(name,value);
     
      goto loop;
    }

  fclose(fp);

  ARGHEAD= ARGLIST;

  /* Now we can use our routines to check & set their own flags ! */

  int input_boolean(char* name,int* value,char* interpret);
  input_boolean("VERBOSE",&i,"off");
  input_boolean("DESCRIBE",&j,"off");
  input_boolean("BEGINNER",&k,"off");
  VERBOSE=i;
  DESCRIBE=j;
  BEGINNER=k;
  
}
#endif

#ifndef __SHUTDOWN_PARSER__
#define __SHUTDOWN_PARSER__
void shutdown_parser(E)
     struct All_variables *E;

{
	if(ARGLIST != NULL) free(ARGLIST);
	if(ARGBUF  != NULL) free(ARGBUF);
	ARGBUF=  NULL;
	ARGLIST= NULL;
	
}
#endif

#ifndef __ADD_TO_PARAMETER_LIST__
#define __ADD_TO_PARAMETER_LIST__
int add_to_parameter_list(name,value)	/* add an entry to arglist, expanding memory */
     register char *name, *value;	/* if necessary */
{
  struct arglist *alptr;
  int len;
  register char *ptr;

  /* check arglist memory */
  if(NLIST >= LISTMAX)
    { LISTMAX += LISTINC;
      if(ARGLIST == NULL)
	ARGLIST= (AL *)malloc(LISTMAX * sizeof(AL));
      else	
	ARGLIST= (AL *)realloc(ARGLIST,LISTMAX * sizeof(AL));
    }
  /* check argbuf memory */
  len= strlen(name) + strlen(value) + 2; /* +2 for terminating nulls */
  if(NBUF+len >= BUFMAX)
    { BUFMAX += BUFINC;
      if(ARGBUF == NULL)
	ARGBUF= (char *)malloc(BUFMAX);
      else	ARGBUF= (char *)realloc(ARGBUF,BUFMAX);
    }
  if(ARGBUF == NULL || ARGLIST == NULL)
   fprintf(stderr,"cannot allocate memory\n");

  /* add name */
  alptr= ARGLIST + NLIST;
  int compute_parameter_hash_table(register char*);
  alptr->hash= compute_parameter_hash_table(name);
  alptr->argname_offset = NBUF;
  ptr= ARGBUF + NBUF;
  do 
    *ptr++ = *name; 
  while(*name++);
  
  /* add value */
  NBUF += len;
  alptr->argval_offset= ptr - ARGBUF;
  do
    *ptr++ = *value;
  while(*value++);
  NLIST++;

}
#endif

#ifndef __COMPUTE_PARAMETER_HASH_TABLES__
#define __COMPUTE_PARAMETER_HASH_TABLES__
compute_parameter_hash_table(s)
     register char *s;
{ register int h;
  
  h= s[0];
  if(s[1])
    h |= (s[1])<<8;
  else
    return(h);
  if(s[2])
    h |= (s[2])<<16;
  else 
    return(h);
  if(s[3])
    h |= (s[3])<<24;
  return(h);
}
#endif

//#ifndef __INPUT_INT__
//#define __INPUT_INT__
int input_int(name,value,interpret)
     char *name;
     int *value;
     char *interpret;

{
    int interpret_control_string();
    struct arglist *alptr; 
    int h, found;
    char  *str;
    
  int exists,essential;
  static double *Default,*minvalue,*maxvalue;
  if(DESCRIBE)
    fprintf(stderr,"input_int: searching for '%s' with default/range '%s'\n",
	    name,(interpret == NULL) ? "**EMPTY**" : interpret);
 
  exists = interpret_control_string(interpret,&essential,&Default,&minvalue,&maxvalue);
  
  if(Default != NULL)
    *value = (int)(*Default);

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;

      str= ARGBUF + alptr->argval_offset;
      sscanf(str,"%d",value);
      found=1;
      break;
    } 

  if(essential && !found)
    { fprintf(stderr,"There MUST be an entry for the parameter %s\n",name);
      exit(12);
    }
  if((minvalue!=NULL) && (*value < (int) *minvalue))
     { *value = (int) *minvalue;
     }
  if((maxvalue!=NULL) && (*value > (int) *maxvalue))
    {  *value = (int) *maxvalue;
    }

  if(VERBOSE)
   { if (found)
       fprintf(stderr,"%25s: (int) = %d \n",name,*value); 
     else
       if (Default != NULL)
	  fprintf(stderr,"%25s: (int) = not found (%d) \n",name,(int)(*Default)); 
       else
	 { fprintf(stderr,"%25s: (int) = not found (no default) \n",name); 
	   if(BEGINNER)
	     { fprintf(stderr,"\t\t Previously set value gives ...");
	       fprintf(stderr,"%d\n",*value);
	     }
	  } 
   }
  //char* substring;
  //substring = strtok(interpret,",");
  if (Default) { free(Default); }
  //substring = strtok(NULL,",");
  if (minvalue) { free(minvalue); }
  //substring = strtok(NULL,",");
  if (maxvalue) { free(maxvalue); }
  return(found);
}
//#endif

#ifndef __INPUT_STRING__
#define __INPUT_STRING__
int input_string(name,value,Default)  /* in the case of a string default=NULL forces input */
     char *name;
     char *value;
     char *Default;
{ 
    char *sptr;
  struct arglist *alptr; 
  int h, hno, hyes, found;
  char line[MAXLINE], *str, *noname;
  int essential;

 
  if(DESCRIBE)
    fprintf(stderr,"input_string: searching for '%s' with default '%s'\n",
	    name,(Default == NULL) ? "no default" : Default);
 
  h=compute_parameter_hash_table(name);
  essential=found=0;

    
    if (Default != NULL)   /* Cannot use "Essential" as this is a valid input */
	strcpy(value,Default);  
    else
	essential=1;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;

      str= ARGBUF + alptr->argval_offset;
      strcpy(value,str);
      found=1;
      break;
    } 
  
  if(essential && !found)
    { fprintf(stderr,"There MUST be an entry for the parameter %s\n",name);
      exit(12);
    }
 
  if(VERBOSE)
    fprintf(stderr,"%25s: (string) = %s (%s)\n",name,
	    (found ? value : "not found"),
	    (Default != NULL ?  Default : "no default")); 

  return(found);
}
#endif

#ifndef __INPUT_BOOLEAN__
#define __INPUT_BOOLEAN__
int input_boolean(name,value,interpret)  /* supports name=on/off too */
     char *name;
     int *value;
     char *interpret;
     
{ char *sptr;
  struct arglist *alptr; 
  int h, hno, hyes, found;
  int interpret_control_string();
  char line[MAXLINE], *str, *noname;

  int essential;
  double *Default,*minvalue,*maxvalue;

  if(DESCRIBE)
    fprintf(stderr,"input_boolean: searching for '%s' with default/range '%s'\n",
	    name,(interpret == NULL) ? "**EMPTY**" : interpret);
 
  interpret_control_string(interpret,&essential,&Default,&minvalue,&maxvalue);

  if(Default != NULL)
    *value = (int)(*Default);
 
  h=compute_parameter_hash_table(name);
  found=0;


  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;

      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 
 
 
  if(!found)
    { if(VERBOSE)
	if (Default != NULL)
	  fprintf(stderr,"%25s: (boolean int) = not found (%d) \n",name,(int)(*Default)); 
	else
	 { fprintf(stderr,"%25s: (boolean int) = not found (no default) \n",name); 
	   if(BEGINNER)
	     { fprintf(stderr,"\t\t Previously set value gives ...");
	       fprintf(stderr,"%d\n",*value);
	     }
	 }
	 
      return(0);
    }
 

  if((strstr(str,"on")!=NULL) || (strstr(str,"ON")!=NULL)) 
    *value=1;
  else if ((strstr(str,"off") != NULL) || (strstr(str,"OFF")!=NULL))
    *value=0;
  else 
    *value=atoi(str);

  if(VERBOSE)
    fprintf(stderr,"%25s: (boolean int) = %d \n",name,*value); 
  
  return(found);
}
#endif

#ifndef __INPUT_FLOAT__
#define __INPUT_FLOAT__
int input_float(name,value,interpret)
     char *name;
     float *value;
     char *interpret;

{ char *sptr;
  struct arglist *alptr;
  
  int h, hno, hyes, found;
  char line[MAXLINE], *str, *noname;
  int exists,essential;
  double *Default,*minvalue,*maxvalue;

 
  if(DESCRIBE)
    fprintf(stderr,"input_float: searching for '%s' with default/range '%s'\n",
	    name,(interpret == NULL) ? "**EMPTY**" : interpret);
 
 
  int interpret_control_string(char*,int*,double**,double**,double**);
  exists=interpret_control_string(interpret,&essential,&Default,&minvalue,&maxvalue);
 
  if(Default != NULL)
    *value = (float) *Default;

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;

      sscanf(str,"%f",value);
      found=1;
      break;
    } 
 
  if(essential && !found)
    { fprintf(stderr,"There MUST be an entry for the parameter %s\n",name);
      exit(12);
    }

  if((minvalue!=NULL) && (*value < (float) *minvalue))
    *value = (float) *minvalue;
  if((maxvalue!=NULL) && (*value > (float) *maxvalue))
    *value = (float) *maxvalue;

  if(VERBOSE)
   { if (found)
       fprintf(stderr,"%25s: (float) = %f \n",name,*value); 
     else
       if (Default != NULL)
	  fprintf(stderr,"%25s: (float) = not found (%f) \n",name,*Default); 
       else
	 { fprintf(stderr,"%25s: (float) = not found (no default) \n",name); 
	   if(BEGINNER)
	     { fprintf(stderr,"\t\t Previously set value gives ...");
	       fprintf(stderr,"%g\n",*value);
	     }
	 }
   }
  return(found);
}
#endif
  
#ifndef __INPUT_DOUBLE__
#define __INPUT_DOUBLE__
int input_double(name,value,interpret)
     char *name;
     double *value;
     char *interpret;

{ char *sptr;
  struct arglist *alptr;
  
  int h, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  int exists,essential;
  double *Default,*minvalue,*maxvalue;


  if(DESCRIBE)
   fprintf(stderr,"input_double: searching for '%s' with default/range '%s'\n",
	   name,(interpret == NULL) ? "**EMPTY**" : interpret);
 
 
  int interpret_control_string(char*,int*,double**,double**,double**);
  exists=interpret_control_string(interpret,&essential,&Default,&minvalue,&maxvalue);
 
  if(Default != NULL)
    *value = *Default;

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      sscanf(str,"%lf",value);
      found=1;
      break;
    } 
 
  if(essential && !found)
    { fprintf(stderr,"There MUST be an entry for the parameter %s\n",name);
      exit(12);
    }
  if((minvalue!=NULL) && (*value <  *minvalue))
    *value =  *minvalue;
  if((maxvalue!=NULL) && (*value >  *maxvalue))
    *value =  *maxvalue;

  if(VERBOSE)
   { if (found)
       fprintf(stderr,"%25s: (double) = %g \n",name,*value); 
     else
       if (Default != NULL)
	  fprintf(stderr,"%25s: (double) = not found (%g) \n",name,*Default); 
       else
	  { fprintf(stderr,"%25s: (double) = not found (no default)\n",name); 
	    if(BEGINNER)
	       { fprintf(stderr,"\t\t Previously set value gives ...");
		 fprintf(stderr,"%g\n",*value);
	       }
	  }
   }
  

  return(found);
}
#endif

#ifndef __INPUT_INT_VECTOR__
#define __INPUT_INT_VECTOR__
int input_int_vector(name,number,value)
     char *name;
     int number;
     int *value; /* comma-separated list of ints */

{ char *sptr;
  struct arglist *alptr;
  char control_string[500];
 
  int h,i, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  if(DESCRIBE)
    fprintf(stderr,"input_int_vector: searching for %s (%d times)\n",name,number);

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 
  /* now interpret vector */
  
  if(!found) return(0);

  for(h=0;h<number;h++)
    { sprintf(control_string,"%s","");
      for(i=0;i<h;i++)
	strcat(control_string,"%*f,");
      strcat(control_string,"%d");
      sscanf(str,control_string,&(value[h]));
    }

  if(VERBOSE)
   fprintf(stderr,"%25s: (vector) = %s\n",name,str); 

  return(found);
}
#endif

#ifndef __INPUT_CHAR_VECTOR__
#define __INPUT_CHAR_VECTOR__
int input_char_vector(name,number,value)
     char *name;
     int number;
     char *value; /* comma-separated list of ints */

{ char *sptr;
  struct arglist *alptr;
  char control_string[500];
 
  int h,i, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  if(DESCRIBE)
    fprintf(stderr,"input_char_vector: searching for %s (%d times)\n",name,number);

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 
  /* now interpret vector */
  
  if(!found) return(0);

  for(h=0;h<number;h++)
    { sprintf(control_string,"%s","");
      for(i=0;i<h;i++)
	strcat(control_string,"%*c,");
      strcat(control_string,"%c");
      sscanf(str,control_string,&(value[h]));
    }

  if(VERBOSE)
   fprintf(stderr,"%25s: (vector) = %s\n",name,str); 

  return(found);
}
#endif

#ifndef __INPUT_FLOAT_VECTOR__
#define __INPUT_FLOAT_VECTOR__
int input_float_vector(char* name, int number, const float* value) {
  char *sptr;
  struct arglist *alptr;
  char control_string[500];
 
  int h,i, hno, hyes, found;
  char line[MAXLINE], *str, *noname;

  if (!number || number == 0) {
    return(0);
  }

  if(DESCRIBE)
    fprintf(stderr,"input_float_vector: searching for %s (%d times)\n",name,number);

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 
  /* now interpret vector */
  
  if(!found) return(0);

  for(h=0;h<number;h++)
    { sprintf(control_string,"%s","");
      for(i=0;i<h;i++)
	strcat(control_string,"%*f,");
      strcat(control_string,"%f");
      sscanf(str,control_string,&(value[h]));
    }

  if(VERBOSE)
   fprintf(stderr,"%25s: (float vector) = %s\n",name,str); 

  return(found);
}
#endif

#ifndef __INPUT_DOUBLE_VECTOR__
#define __INPUT_DOUBLE_VECTOR__
int input_double_vector(name,number,value)
     char *name;
     int number;
     double *value; /* comma-separated list of floats */

{ char *sptr;
  struct arglist *alptr;
  char control_string[500];
 
  int h,i, hno, hyes, found;
  char line[MAXLINE], *str, *noname;
 
  if(DESCRIBE)
    fprintf(stderr,"input_double_vector: searching for %s (%d times)\n",name,number);

  h=compute_parameter_hash_table(name);
  found=0;

  /* search list backwards, stopping at first find */
  for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
    { if(alptr->hash != h)
	continue;
      if(strcmp(ARGBUF+alptr->argname_offset,name) != 0)
	continue;
      str= ARGBUF + alptr->argval_offset;
      found=1;
      break;
    } 

  if(!found) return(0);

 /* now interpret vector */
  
  for(h=0;h<number;h++)
    { sprintf(control_string,"%s","");
      for(i=0;i<h;i++)
	strcat(control_string,"%*f,");
      strcat(control_string,"%lf");
      sscanf(str,control_string,&(value[h]));
    }

  if(VERBOSE)
   fprintf(stderr,"%25s: (double vector) = %s\n",name,str); 

  return(found);
}
#endif

/* =================================================== */

#ifndef __INTERPRET_CONTROL_STRING__
#define __INTERPRET_CONTROL_STRING__
int interpret_control_string(interpret,essential,Default,minvalue,maxvalue)
     char *interpret;
     int *essential;
     double **Default,**minvalue,**maxvalue;

{ char *substring;

  *Default=*maxvalue=*minvalue=NULL;
  *essential=0;
     
    return(0);  /* nothing to interpret */

  if ((substring=strtok(interpret,",")) == NULL)
    return(0);  /* nothing to interpret */
    
 
  if (strstr(substring,"essential")!=NULL)
    *essential=1; /* no default possible, must read a value */
  else
    if (strstr(substring,"nodefault")==NULL) 
     { *Default = (double *) malloc(sizeof(double));
       if((strstr(substring,"on")!=NULL) || (strstr(substring,"ON")!=NULL))
	 **Default = 1.0;
       else 
	 if ((strstr(substring,"off") != NULL) || (strstr(substring,"OFF")!=NULL))
	   **Default = 0.0; 
	 else
	   sscanf(substring,"%lf",*Default);  /* read number as a default value */
       
     }
  
  if ((substring=strtok(NULL,",")) == NULL) /* minvalue */
    { /* no minimum, no maximum */
      return(1);
    }

  if (strstr(substring,"nomin")==NULL)
    { *minvalue = (double *) malloc(sizeof(double));
      sscanf(substring,"%lf",*minvalue);
    }
  
  if ((substring=strtok(NULL,",")) == NULL) /* maxvalue */
    { /* no maximum */
      if (DESCRIBE)
	fprintf(stderr,"minimum but no maximum\n");
      return(1);
    }

  if (strstr(substring,"nomax")==NULL)
    { *maxvalue = (double *) malloc(sizeof(double));
      sscanf(substring,"%lf",*maxvalue);
    }

  return(1);
 
}
#endif
