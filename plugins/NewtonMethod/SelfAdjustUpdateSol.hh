// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_SelfAdjustUpdateSol_hh
#define COOLFluiD_Numerics_NewtonMethod_SelfAdjustUpdateSol_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonMethod/StdUpdateSol.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class SelfAdjustUpdateSol : public StdUpdateSol {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit SelfAdjustUpdateSol(const std::string& name);

  /// Destructor.
  virtual ~SelfAdjustUpdateSol()
  {
  }

  /// Execute Processing actions
  virtual void execute();

protected:

}; // class SelfAdjustUpdateSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_SelfAdjustUpdateSol_hh
