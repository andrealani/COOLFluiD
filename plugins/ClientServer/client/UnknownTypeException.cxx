#include <string>

#include "Common/CodeLocation.hh"

#include "ClientServer/client/UnknownTypeException.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::client;

UnknownTypeException::UnknownTypeException(const CodeLocation& where, 
                                           const std::string& what) 
: Exception(where, what, "UnknownTypeException") 
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

UnknownTypeException::UnknownTypeException(const UnknownTypeException& e) throw() 
: Exception(e) 
{
  
}