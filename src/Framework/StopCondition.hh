// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_StopCondition_hh
#define COOLFluiD_Framework_StopCondition_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
     class ConvergenceStatus;

//////////////////////////////////////////////////////////////////////////////

/// This class defines the condition to be satisfied to stop
/// the computation
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API StopCondition :
    public Common::OwnedObject,
    public Config::ConfigObject,
    public Common::NonCopyable<StopCondition> {

public:

  typedef Environment::ConcreteProvider<StopCondition,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Constructor
  StopCondition(const std::string& name);

  /// Default destructor
  virtual ~StopCondition();

  /// Gets the Class name
  static std::string getClassName() { return "StopCondition"; }

  /// returns true if we have to evaluate the stopcondition on all CPU domains.
  virtual bool IsGlobal () const = 0;

  /// Returns true if the SubSystem should end
  virtual bool isAchieved (const ConvergenceStatus& status) = 0;

}; // end of class StopCondition

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(StopCondition) // define the factory

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_StopCondition_hh
