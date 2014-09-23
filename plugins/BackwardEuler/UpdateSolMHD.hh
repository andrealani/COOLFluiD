// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_BackwardEuler_UpdateSolMHD_hh
#define COOLFluiD_Numerics_BackwardEuler_UpdateSolMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "BwdEulerData.hh"
#include "Framework/DataSocketSink.hh"
#include "MHD/MHDTerm.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Physics::MHD;

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class UpdateSolMHD : public BwdEulerCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit UpdateSolMHD(const std::string& name);

  /// Destructor.
  ~UpdateSolMHD()
  {
  }

  /// Set up the private data
  void setup();

  /// Execute Processing actions
  void execute();

  /// Configures this object by complementing the
  /// implementation in ConfigObject
  void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /// socket for Rhs
  Framework::DataSocketSink<
                              CFreal> socket_rhs;

  /// safe pointer to the physical data
  Common::SafePtr<MHDTerm> _model;

  /// pressure correction value
  CFreal _pressureCorrectionVal;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<
                              CFreal> socket_updateCoeff;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

}; // class UpdateSolMHD

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_BackwardEuler_UpdateSolMHD_hh

