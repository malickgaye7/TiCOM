#ifndef __PARSING_H__
#define __PARSING_H__

/* A new header file created to resolve symbol duplication bugs.
   Last updated 1/17/2022 */

#include <stdio.h>

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

//#include <malloc.h>
#include <sys/types.h>
#ifndef __sunos__
#include <strings.h>
#endif

#ifndef STRING_H
#define STRING_H
#include <string.h>
#endif

#ifndef __GLOBAL_DEFS_H__
#define __GLOBAL_DEFS_H__
#include "global_defs.h"
#endif

#define MAXLINE		1024	/* max length of line in input file */
#define MAXNAME		64	/* max length of name */
#define MAXVALUE	1024	/* max length of value */
#define MAXFILENAME	64	/* max length of par file name */
#define MAXVECTOR	10	/* max # of elements for unspecified vectors */

/* abbreviations: */
#define AL 		struct arglist
#define PROGNAME	ext_par.progname
#define FLAGS		ext_par.argflags
#define ARGLIST		ext_par.arglist
#define ARGHEAD		ext_par.arghead
#define ARGBUF		ext_par.argbuf
#define NLIST		ext_par.nlist
#define NBUF		ext_par.nbuf
#define LISTMAX		ext_par.listmax
#define BUFMAX		ext_par.bufmax
#define LISTFILE	ext_par.listout

#define LISTINC		32	/* increment size for arglist */
#define BUFINC		1024	/* increment size for argbuf */

extern int VERBOSE;
extern int DESCRIBE;
extern int BEGINNER;

//extern struct ext_par;
//extern struct arglist;

void setup_parser(struct All_variables*, int, char**);
void shutdown_parser(struct All_variables*);
int add_to_parameter_list(register char*, register char*);
int compute_parameter_hash_table(register char*);
int input_int(char*, int*, char*);
int input_string(char*, char*, char*);
int input_boolean(char*, int*, char*);
int input_float(char*, float*, char*);
int input_double(char*, double*, char*);
int input_int_vector(char*, int, int*);
int input_char_vector(char*, int, char*);
int input_float_vector(char*, int, const float*);
int input_double_vector(char*, int, double*);
int interpret_control_string(char*, int*, double**, double**, double**);

#endif /* __PARSING_H__ */