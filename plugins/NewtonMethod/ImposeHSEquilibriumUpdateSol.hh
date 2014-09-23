// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_ImposeHSEquilibriumUpdateSol_hh
#define COOLFluiD_Numerics_NewtonMethod_ImposeHSEquilibriumUpdateSol_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdUpdateSol.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class ImposeHSEquilibriumUpdateSol : public StdUpdateSol {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit ImposeHSEquilibriumUpdateSol(const std::string& name);

  /// Destructor.
  virtual ~ImposeHSEquilibriumUpdateSol()
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
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
protected:

  /// storage of volumes
  Framework::DataSocketSink < CFreal > socket_volumes;
  
  /// storage of gravity 
  Framework::DataSocketSource < CFreal > socket_gravity;
    
}; // class ImposeHSEquilibriumUpdateSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_ImposeHSEquilibriumUpdateSol_hh
