#include <string>

#include "Common/CodeLocation.hh"

#include "ClientServer/server/UnknownClientIdException.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::server;

UnknownClientIdException::UnknownClientIdException(const CodeLocation& where,
                                                   const std::string& what) 
: Exception(where, what, "UnknownClientIdException") 
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

UnknownClientIdException::UnknownClientIdException(const UnknownClientIdException& e) throw() 
: Exception(e) 
{
  
}