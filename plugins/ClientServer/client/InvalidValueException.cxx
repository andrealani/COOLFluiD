#include <string>

#include "Common/CodeLocation.hh"

#include "ClientServer/client/InvalidValueException.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::client;

InvalidValueException::InvalidValueException(const CodeLocation& where, 
                                             const std::string& what) 
: Exception(where, what, "InvalidValueException") 
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

InvalidValueException::InvalidValueException(const InvalidValueException& e) throw() 
: Exception(e) 
{
  
}