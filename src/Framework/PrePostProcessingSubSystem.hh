// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PrePostProcessingSubSystem_hh
#define COOLFluiD_Framework_PrePostProcessingSubSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StandardSubSystem.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// A PrePostProcessingSubSystem is a concrete implementation of the SubSystem
/// that only apply a pre/post processing.
/// @author Andrea Lani
class Framework_API PrePostProcessingSubSystem : public StandardSubSystem {

public: // functions

  /// The default constructor without arguments.
  /// @see ~PrePostProcessingSubSystem()
  PrePostProcessingSubSystem(const std::string& name);

  /// Destructor
  /// @see PrePostProcessingSubSystem()
  ~PrePostProcessingSubSystem();

  /// Setup Phase. Setup parameters of the PrePostProcessingSubSystem.
  /// @see MeshData
  virtual void setup();

  /// Run (Process) Phase. All the big number crunching work goes here.
  virtual void run();

  /// Unsetup Phase.
  virtual void unsetup();

  /// Configures this Simualtion.
  /// Sets up the data for this Object.
  virtual void configure( Config::ConfigArgs& args );

protected: // functions

  /// Gets a vector with all the Methods this PrePostProcessingSubSystem has
  /// @return vector with Method pointers.
  std::vector<Framework::Method*> getMethodList();
  
  /// Set the Method's collaborators.
  void setCollaborators();

}; // class PrePostProcessingSubSystem

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PrePostProcessingSubSystem_hh
