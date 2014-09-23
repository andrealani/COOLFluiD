#include <string>

#include "Common/CodeLocation.hh"

#include "ClientServer/server/TypesNotFoundException.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::server;

TypesNotFoundException::TypesNotFoundException(const CodeLocation& where,
                                                   const std::string& what) 
 : Exception(where, what, "UnknownClientIdException") 
{

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TypesNotFoundException::TypesNotFoundException(const TypesNotFoundException& e) throw() 
 : Exception(e) 
{

}