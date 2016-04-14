// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_UpdateSolFVMCC_hh
#define COOLFluiD_Numerics_NewtonMethod_UpdateSolFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonMethod/StdUpdateSol.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class UpdateSolFVMCC : public StdUpdateSol {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit UpdateSolFVMCC(const std::string& name);

  /// Destructor.
  virtual ~UpdateSolFVMCC()
  {
  }

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:
  
  /// Correct the given unphysical state
  virtual void correctUnphysicalStates(const std::vector<CFuint>& badStatesIDs);
  
protected:
  
  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  /// corrected state
  RealVector m_correctedState;

}; // class UpdateSolFVMCC

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_UpdateSolFVMCC_hh
