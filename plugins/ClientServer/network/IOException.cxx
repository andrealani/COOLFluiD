#include <string>

#include "Common/CodeLocation.hh"

#include "ClientServer/network/IOException.hh"

using namespace COOLFluiD::Common;
using namespace COOLFluiD::network;

IOException::IOException(const CodeLocation& where, const std::string& what) 
: Exception(where, what, "IOException") 
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IOException::IOException(const IOException& e) throw() 
: Exception(e) 
{
  
}