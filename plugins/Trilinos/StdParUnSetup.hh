// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_StdParUnSetup_hh
#define COOLFluiD_Numerics_Trilinos_StdParUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Trilinos.hh"
#include "TrilinosLSSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to deallocate data specific to the
 * method to which it belongs.
 *
 * @author Tim Boonen
 */
class StdParUnSetup : public TrilinosLSSCom {
public:

  /**
   * Constructor.
   */
  explicit StdParUnSetup(const std::string& name) : TrilinosLSSCom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdParUnSetup()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

}; // class UnSeqSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Trilinos_StdParUnSetup_hh

