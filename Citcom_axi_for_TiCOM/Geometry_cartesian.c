#include "element_definitions.h"
#include "global_defs.h"



void set_2dc_defaults(E)
     struct All_variables *E;
{ 

  E->mesh.nsd = 2;
  E->mesh.dof = 2;
  
}


void set_2pt5dc_defaults(E)  
    struct All_variables *E;
{ 

  E->mesh.nsd = 2;
  E->mesh.dof = 3;
 
}

void set_3dc_defaults(E)
     struct All_variables *E;
{ 

  E->mesh.nsd = 3;
  E->mesh.dof = 3;
 
}
