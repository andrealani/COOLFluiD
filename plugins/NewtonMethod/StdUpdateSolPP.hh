// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_StdUpdateSolPP_hh
#define COOLFluiD_Numerics_NewtonMethod_StdUpdateSolPP_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonMethod/NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class StdUpdateSolPP : public NewtonIteratorCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit StdUpdateSolPP(const std::string& name);

  /// Destructor.
  virtual ~StdUpdateSolPP()
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
  
  /// Returns the DataSocket's that this command needs as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
protected:
  
  /// Correct the given unphysical state
  virtual void correctUnphysicalStates(const std::vector<CFuint>& badStatesIDs) {} 
  
protected:

  /// handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// handle to rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// handle to update coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// handle to invalid states
  Framework::DataSocketSource<CFreal> socket_invalidStates;

  /// the relaxation parameter vector
  std::vector<CFreal> m_alpha;
  
  /// Check that each update creates variables with physical meaning
  bool m_validate;

  /// boundary pressure value
  CFreal _pBC;
  /// boundary plasma density value
  CFreal _rhoBC;

}; // class StdUpdateSolPP

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_StdUpdateSolPP_hh
