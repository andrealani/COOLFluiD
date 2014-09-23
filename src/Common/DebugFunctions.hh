// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_DebugFunctions_hh
#define COOLFluiD_Common_DebugFunctions_hh

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// pause and prompt for input user till the specified character is given
template <int NUM>
static void DEBUG_CONDITIONAL_PAUSE(char stop = '$')
{
  static char in; 
  if (PE::GetPE().GetRank() == 0) {
    std::cout << "DEBUG_CONDITIONAL_PAUSE: input a character different from \"" << stop << "\" to stop again \n";
    if (in != stop) {std::cin >> in;}
  }
  
#ifdef CF_HAVE_MPI
  char tmp = in; 
  MPI_Bcast(&tmp, 1, MPI_CHAR, 0, PE::GetPE().GetCommunicator());
  in = tmp;
#endif
}
      
//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_DebugFunctions_hh
