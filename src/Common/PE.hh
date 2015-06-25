// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_PE_hh
#define COOLFluiD_Common_PE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/PEInterface.hh"
#include "Common/WorkerStatus.hh"

#ifdef CF_HAVE_MPI
#include "Common/MPI/PEInterfaceMPI.hh"
#include "Common/MPI/MPIError.hh"
#else
#include "Common/SERIAL/PEInterfaceSERIAL.hh"
#endif // CF_HAVE_MPI

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class controls the Common environment
/// @author Dries Kimpe
class Common_API PE {
public:
  
  /// Initialise the PE
  static void InitPE (int* argc, char*** args);
  
  /// Checks if the PE is initialized
  static bool IsInitialised ();
  
  /// Return a reference to the current PE
  static PEInterface<>& GetPE ();
  
  /// Free the PE
  static void DonePE ();
  
 private:
  
  /// the current PE
  static PEInterface<>* m_curr_PE;
  
  /// Flag to keep track of Common Enviroment Initialized
  /// Note: cannot rely on m_curr_PE pointer because objects held by the enviroment
  ///       could not then check for initialization. An egg and chicken problem,
  ///       that appeared when using CFLog on the destructors of PEInterface related objects.
  static bool m_is_init;
       
}; // end class PE

//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD

#endif // COOLFluiD_Common_PE_hh
