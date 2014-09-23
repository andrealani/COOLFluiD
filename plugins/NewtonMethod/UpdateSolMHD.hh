// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_UpdateSolMHD_hh
#define COOLFluiD_Numerics_NewtonMethod_UpdateSolMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"
#include "MHD/MHDTerm.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Physics::MHD;

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class UpdateSolMHD : public NewtonIteratorCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit UpdateSolMHD(const std::string& name);

  /// Destructor.
  virtual ~UpdateSolMHD()
  {
  }

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Execute Processing actions
  virtual void execute();

  /// Configures this object by complementing the
  /// implementation in ConfigObject
  void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// handle to update coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// safe pointer to the physical data
  Common::SafePtr<MHDTerm> _model;

  /// pressure correction value
  CFreal _pressureCorrectionVal;

  /// the relaxation parameter vector
  std::vector<CFreal> m_alpha;
  
  /// Check that each update creates variables with physical meaning
  bool m_validate;

}; // class UpdateSolMHD

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_UpdateSolMHD_hh
