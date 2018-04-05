// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeTerm_hh
#define COOLFluiD_Framework_ComputeTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Environment/ConcreteProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Common/NotImplementedException.hh"
#include "MethodStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class VarSet;
    class GeometricEntity;

//////////////////////////////////////////////////////////////////////////////

/// This class offers a basic interface for all flux splitters
/// @author Andrea Lani
/// @author Tiago Quintino

template < typename METHODDATA>
class ComputeTerm : public MethodStrategy<METHODDATA> {
public:

  typedef Environment::ConcreteProvider<ComputeTerm> PROVIDER;

  /// Constructor
  ComputeTerm(const std::string& name) :
    MethodStrategy<METHODDATA>(name)
  {
  }
  
  /// Default destructor
  virtual ~ComputeTerm()
  {
  }

  /// Set up private data to prepare the simulation
  virtual void setup() 
  {
    MethodStrategy<METHODDATA>::setup();
  }
  
  /// Compute the term in the current cell
  virtual void computeTerm(Framework::GeometricEntity* const cell, RealMatrix& result)
  {
    throw Common::NotImplementedException (FromHere(),getClassName() + "::computeTerm(GeometricEntity*, RealMatrix)");
  }

  /// Compute the term in the current cell
  virtual void computeTerm(Framework::GeometricEntity* const cell, RealVector& result)
  {
    throw Common::NotImplementedException (FromHere(),getClassName() + "::computeTerm(GeometricEntity*, RealVector)");
  }

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Config::ConfigObject::configure(args);
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeTerm";
  }

private: // methods

  /// Private Copy Constructor
  ComputeTerm(const ComputeTerm& f);

  /// Private Assignement operator
  const ComputeTerm& operator=(const ComputeTerm& f);

protected: // data

  /// i state
  CFuint _iState;

  /// j state
  CFuint _jState;

}; // end of class ComputeTerm

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeTerm_hh
