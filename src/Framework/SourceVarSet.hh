// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SourceVarSet_hh
#define COOLFluiD_Framework_SourceVarSet_hh

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

/// This class represents the basic interface for a Source
/// variable set, which provides physical model dependent data
/// and methods associated to a choice of variables
/// @author Andrea Lani
class Framework_API SourceVarSet : public Common::OwnedObject,
                     public Config::ConfigObject,
                     public Common::NullableObject {
public:

  typedef Environment::ConcreteProvider<SourceVarSet,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Constructor
  SourceVarSet(const std::string& name);

  /// Default destructor
  virtual ~SourceVarSet();

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

  /// Checks if the source term has a part independent from the solution
  virtual bool hasIndepCoef() const = 0;

  /// Checks if the source term has a linear part. By default it doesnt.
  virtual bool hasLinearCoef() const = 0;

  /// Gets the Independent Source Coeficients
  virtual void getIndepSourceCoefs(const Framework::State& state, const RealVector& normals, RealVector& coef)
  {
    throw Common::NotImplementedException (FromHere(),getClassName());
  }

  /// Gets the Independent Source Coeficients
  virtual void getIndepSourceCoefs(const Framework::State& state, RealVector& coef)
  {
    throw Common::NotImplementedException (FromHere(),getClassName());
  }

  /// Gets the Linear Source Coeficients
  virtual void getLinearSourceCoefs(const Framework::State& state, const RealVector& normals, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),getClassName());
  }

  /// Gets the Linear Source Coeficients
  virtual void getLinearSourceCoefs(const Framework::State& state, RealMatrix& coef)
  {
    throw Common::NotImplementedException (FromHere(),getClassName());
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "SourceVarSet";
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
  SourceVarSet(const SourceVarSet& v);

  /// Private Assignement operator
  const SourceVarSet& operator=(const SourceVarSet& v);

protected:

  /// list of the variables names
  std::vector<std::string> _varNames;

}; // end of class SourceVarSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(SourceVarSet)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SourceVarSet_hh
