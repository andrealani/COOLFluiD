// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeFlux_hh
#define COOLFluiD_Framework_ComputeFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/MethodStrategy.hh"
#include "MathTools/RealVector.hh"
#include "Environment/ConcreteProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class ConvectiveVarSet;
    class GeometricEntity;

//////////////////////////////////////////////////////////////////////////////

/// This class offers a basic interface for all flux splitters
/// @author Andrea Lani
/// @author Tiago Quintino
template < typename METHODDATA >
class ComputeFlux : public Framework::MethodStrategy<METHODDATA> {
public:

  typedef Environment::ConcreteProvider<ComputeFlux,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Constructor
  ComputeFlux(const std::string& name) :
    MethodStrategy<METHODDATA>(name)
  {
  }
  
  /// Default destructor
  virtual ~ComputeFlux()
  {
  }
  
  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Config::ConfigObject::configure(args);
  }
  
  /// Set up private data to prepare the simulation
  virtual void setup() 
  {
    Framework::MethodStrategy<METHODDATA>::setup();
  }
  
  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<BaseDataSocketSink> > result;
    return result;
  }
  
  /// Compute the flux in the current face
  virtual void computeFlux(RealVector& result) = 0;
  
  /// Gets the Class name
  static std::string getClassName() {  return "ComputeFlux"; }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
private: // methods

  /// Private Copy Constructor
  ComputeFlux(const ComputeFlux& f);

  /// Private Assignement operator
  const ComputeFlux& operator=(const ComputeFlux& f);
  
}; // end of class ComputeFlux

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeFlux_hh
