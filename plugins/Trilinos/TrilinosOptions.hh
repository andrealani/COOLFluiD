// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_TrilinosOptions_hh
#define COOLFluiD_Numerics_Trilinos_TrilinosOptions_hh

//////////////////////////////////////////////////////////////////////////////


#include "AztecOO.h"

#include "Trilinos.hh"
#include "Common/StringOps.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

/**
 *  Class for mapping descriptions of options of the Trilinos solver interface
 *  to their internal codings.
 *
 *  @author   Tim Boonen
 *  @author   Tiago Quintino
 */
class TrilinosOptions {
public:

  /**
   * Set all the options
   */
  static void setAllOptions();

private:
  /**
   * Default constructor without arguments.
   */
  TrilinosOptions();

  /**
   * Destructor.
   */
  ~TrilinosOptions();

 private: // data

}; // end of class TrilinosOptions

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Trilinos_TrilinosOptions_hh
