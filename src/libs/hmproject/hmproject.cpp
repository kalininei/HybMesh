#include "hmproject.h"
#include "Gmsh.h"
#include <libxml/xmlmemory.h>

namespace{

//gmesh library initialization
struct gmsh_initializer{
	gmsh_initializer(){GmshInitialize();}
	~gmsh_initializer(){GmshFinalize();}
} gi;

//This object is used to invoke xmlCleanupParser() at the end of the library utilization
//in order to prevent memory leaks.
struct xmlCleanUp{
	~xmlCleanUp(){ xmlCleanupParser(); }
} xmlcleanup;

}//namespace
