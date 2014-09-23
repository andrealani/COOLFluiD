// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_PEFunctions_hh
#define COOLFluiD_Common_PEFunctions_hh

#include "Common/PE.hh"  

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////
      
/// run serial a pointer to member function taking no arguments
/// @author Andrea Lani
template<typename R, typename C, R(C::*M)()>
static void runSerial(C* obj) 
{
  if (PE::GetPE().IsParallel()) {
    PE::GetPE().setBarrier();
    for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(); ++i) {
      if (i == PE::GetPE().GetRank ()) {
	(obj->*M)();
      }
      PE::GetPE().setBarrier();
    }
  }
  else {
    (obj->*M)();
  }
}
 
//////////////////////////////////////////////////////////////////////////////
      
/// run serial a pointer to member function with arity = 1
/// @author Andrea Lani
template<typename R, typename A1, typename C, R(C::*M)(A1)>
void runSerial(C* obj, A1 arg1) 
{
  if (PE::GetPE().IsParallel()) {
    PE::GetPE().setBarrier();
    
    for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(); ++i) {
      if (i == PE::GetPE().GetRank ()) {
	(obj->*M)(arg1);
      }
      PE::GetPE().setBarrier();
    }
  }
  else {
    (obj->*M)(arg1);
  }
}
 
//////////////////////////////////////////////////////////////////////////////
       
/// run serial a pointer to member function with arity = 2
/// @author Andrea Lani
template<typename R, typename A1, typename A2, typename C, R(C::*M)(A1, A2)>
void runSerial(C* obj, A1 arg1, A2 arg2) 
{
  if (PE::GetPE().IsParallel()) {
    PE::GetPE().setBarrier();
    
    for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(); ++i) {
      if (i == PE::GetPE().GetRank ()) {
	(obj->*M)(arg1, arg2);
      }
      PE::GetPE().setBarrier();
    }
  }
  else {
    (obj->*M)(arg1, arg2);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_PEFunctions_hh
