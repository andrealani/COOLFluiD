// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_UpdateSol_hh
#define COOLFluiD_Numerics_ForwardEuler_UpdateSol_hh

//////////////////////////////////////////////////////////////////////////////

#include "FwdEulerData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class  ForwardEuler_API UpdateSol : public FwdEulerCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit UpdateSol(const std::string& name);
  
  /// Destructor.
  ~UpdateSol()
  {
  }

  /// Execute Processing actions
  virtual void execute();
  
  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:
  
  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// handle to update coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// handle to volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// array of flags for active states (states to be updated)
  std::vector<bool> m_activeStates;

  /// Check that eac h update creates variables with physical meaning
  bool m_validate;
  
  /// Flag to clip residual to 0 if < 1e-16
  bool m_clipResidual;
  
  /// Names of the TRSs whose states have a RHS to be set to 0
  std::vector<std::string> m_trsWithNullRHS;
  
}; // class UpdateSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ForwardEuler_UpdateSol_hh

