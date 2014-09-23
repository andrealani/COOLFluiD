// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_InertiaVarSet_hh
#define COOLFluiD_Framework_InertiaVarSet_hh

//////////////////////////////////////////////////////////////////////////////



#include "Common/NullableObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/State.hh"
#include "Common/OwnedObject.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a Inertia
/// variable set, which provides physical model dependent data
/// and methods associated to a choice of variables
/// @author Andrea Lani
class Framework_API InertiaVarSet : public Common::OwnedObject,
                      public Config::ConfigObject,
                      public Common::NullableObject {
public:

  typedef Environment::ConcreteProvider<InertiaVarSet,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Constructor
  InertiaVarSet(const std::string& name);

  /// Default destructor
  virtual ~InertiaVarSet();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    ConfigObject::configure(args);

    // add here configuration, specific of this class
  }

  /// Setup the member data
  virtual void setup()
  {
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "InertiaVarSet";
  }

  /// Get variable names
  const std::vector<std::string>& getVarNames() const
  {
    cf_assert(_varNames.size() > 0);
    return _varNames;
  }

protected: // methods

  /// Set the list of the variable names
  void setVarNames(const std::vector<std::string>& varNames)
  {
    if (_varNames.size() == 0) {
      _varNames = varNames;
    }
  }

private: // methods

  /// Private Copy Constructor
  InertiaVarSet(const InertiaVarSet& v);

  /// Private Assignement operator
  const InertiaVarSet& operator=(const InertiaVarSet& v);

protected:

  /// list of the variables names
  std::vector<std::string> _varNames;

}; // end of class InertiaVarSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(InertiaVarSet) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_InertiaVarSet_hh
