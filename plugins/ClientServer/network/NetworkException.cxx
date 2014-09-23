#include <string>

#include "Common/CodeLocation.hh"

#include "ClientServer/network/NetworkException.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::network;

NetworkException::NetworkException(const CodeLocation& where, 
                                   const std::string& what) 
: Exception(where, what, "NetworkException") 
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NetworkException::NetworkException(const NetworkException& e) throw() 
: Exception(e) 
{
  
}