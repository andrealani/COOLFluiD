#include <string>

#include "Common/CodeLocation.hh"

#include "ClientServer/network/MalformedDataException.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::network;

MalformedDataException::MalformedDataException(const CodeLocation& where, 
                                               const std::string& what) 
: Exception(where, what, "MalformedDataException") 
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MalformedDataException::MalformedDataException(const MalformedDataException& e) throw() 
: Exception(e) 
{
  
}