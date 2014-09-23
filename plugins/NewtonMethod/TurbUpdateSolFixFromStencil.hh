// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_TurbUpdateSolFixFromStencil_hh
#define COOLFluiD_Numerics_NewtonMethod_TurbUpdateSolFixFromStencil_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class TurbUpdateSolFixFromStencil : public NewtonIteratorCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit TurbUpdateSolFixFromStencil(const std::string& name);

  /// Destructor.
  virtual ~TurbUpdateSolFixFromStencil()
  {
  }

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Execute Processing actions
  virtual void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for stencil
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  // handle to rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // handle to update coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// the relaxation parameter vector
  std::vector<CFreal> _alpha;

  //free stream values of the turbulent variables
  CFreal _KInlet;
  CFreal _OmegaInlet;

}; // class TurbUpdateSolFixFromStencil

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_TurbUpdateSolFixFromStencil_hh
