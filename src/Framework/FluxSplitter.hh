// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FluxSplitter_hh
#define COOLFluiD_Framework_FluxSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeFlux.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class ConvectiveVarSet;
    class DiffusiveVarSet;
    class GeometricEntity;

//////////////////////////////////////////////////////////////////////////////

/// This class offers a basic interface for all flux splitters
/// @author Andrea Lani
/// @author Tiago Quintino
template < typename METHODDATA >
class FluxSplitter : public ComputeFlux<METHODDATA> {
public: // typedef

  /// typedef needed for the registration into the factory
  typedef Framework::BaseMethodStrategyProvider<METHODDATA,FluxSplitter<METHODDATA> > PROVIDER;

public: // methods

  /// Constructor
  FluxSplitter(const std::string& name) : ComputeFlux<METHODDATA>(name)
  {
  }

  /// Default destructor
  virtual ~FluxSplitter()
  {
  }

  /// Set up private data to prepare the simulation
  virtual void setup()
  {
    ComputeFlux<METHODDATA>::setup();
    _lFluxJacobian.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
    _rFluxJacobian.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  }
  
  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<BaseDataSocketSink> > result = ComputeFlux<METHODDATA>::needsSockets();
    return result;
  }
  
  /// Compute the flux in the current face
  virtual void computeFlux(RealVector& result) = 0;

  /// Get the flux jacobian of the right state
  virtual Common::SafePtr<RealMatrix> getRightFluxJacob()
  {
    return &_rFluxJacobian;
  }
  
  /// Get the flux jacobian of the left state
  virtual Common::SafePtr<RealMatrix> getLeftFluxJacob()
  {
    return &_lFluxJacobian;
  }
  
  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Config::ConfigObject::configure(args);
  }
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "FluxSplitter";
  }
  
protected:

  /// jacobian matrix of the convective fluxes for the left state
  RealMatrix _lFluxJacobian;

  /// jacobian matrix of the convective fluxes for the right state
  RealMatrix _rFluxJacobian;


}; // end of class FluxSplitter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FluxSplitter_hh
