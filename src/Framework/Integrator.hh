// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_Integrator_hh
#define COOLFluiD_Framework_Integrator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NoSuchValueException.hh"
#include "Common/SelfRegistPtr.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/IntegratorPattern.hh"
#include "Framework/CFQuadrature.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {
    class FactoryRegistry;
  }
  
  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

/// This class provides the base class of an Integrator.
/// Each type of integrator (ContourIntegrator or VolumeIntegrator)
/// has a list of implementor's (IntegratorImpl), one for each
/// GeometricEntity shape.
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename INTEGRATORIMPL>
class Integrator {
private:

  /// ActiveDataType to store INTEGRATORIMPL's
  typedef typename std::vector<Common::SelfRegistPtr<INTEGRATORIMPL> > ActiveDataType;

  /// InnerDataType to store INTEGRATORIMPL's
  typedef typename std::vector<Common::SelfRegistPtr<INTEGRATORIMPL> > InnerDataType;

  /// DataType to store INTEGRATORIMPL's
  typedef std::vector<InnerDataType> PoolDataType;

public:

  /// Constructor
  Integrator();

  /// Destructor
  ~Integrator();

  /// Set up the integrator for the simulation
  void setup();

  /// Unset up the integrator
  void unsetup();

  /// @return the maximum pattern use in all integrations
  IntegratorPattern getMaxIntegratorPattern() const;

  /// Gets the correct solution integrator for the supplied GeometricEntity
  /// according with the configuration of the SubSystem.
  /// @todo First integrator is always chosen. Fix this to be configuration dependent.
  INTEGRATORIMPL* getSolutionIntegrator(const GeometricEntity* const geo) const
  {
    cf_assert(geo->getSolInterpolatorID() < m_activeIntegrators.size());
    cf_assert(m_activeIntegrators[geo->getSolInterpolatorID()].isNotNull());
    return m_activeIntegrators[geo->getSolInterpolatorID()].getPtr();
  }

  /// Gets the correct geometry integrator for the supplied GeometricEntity
  /// according with the configuration of the SubSystem.
  /// @todo First integrator is always chosen. Fix this to be configuration dependent.
  INTEGRATORIMPL* getGeometryIntegrator(const GeometricEntity* const geo) const
  {
    cf_assert(geo->getGeomInterpolatorID() < m_activeIntegrators.size());
    cf_assert(m_activeIntegrators[geo->getGeomInterpolatorID()].isNotNull());
    return m_activeIntegrators[geo->getGeomInterpolatorID()].getPtr();
  }

  /// Sets the quadrature type and integration order for all GeometricEntity's
  void setIntegrationForAllGeo(const CFQuadrature::Type& quadType,
                               const CFPolyOrder::Type& order);

  /// Sets the quadrature type and integration order for GeometricEntity's of CFGeomShape = Hexa
  void setIntegrationForShape(const CFGeoShape::Type shape,
                              const CFQuadrature::Type& quadType,
                              const CFPolyOrder::Type& order);

  /// Set the factory registry
  virtual void setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr) {m_fr = fr;}
  
private: // helper methods

  /// Activates the integrators from the pool
  void activateIntegrators();

  /// Copy constructor
  Integrator(const Integrator&);

  /// assignment operator
  Integrator& operator= (const Integrator&);

protected: // data

  /// check for integrator setup
  bool _setup;

  /// pool of implementors, one list of possible integrators for
  /// each shape of GeometricEntity
  ///  with a certain type of Interpolator
  PoolDataType m_pool;

  /// list of implementors, one for each shape of GeometricEntity
  ///  with a certain type of Interpolator
  ActiveDataType m_activeIntegrators;

  /// map from CFGeoShape::Type to the pair of CFQuadrature::Type and CFPolyOrder::Type
  std::map<CFGeoShape::Type,std::pair<CFQuadrature::Type,CFPolyOrder::Type> > m_integType;
  
  /// factory registry to allow polymorphic creation of objects
  Common::SafePtr<Common::FactoryRegistry> m_fr;
  
}; // end class Integrator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Integrator.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_Integrator_hh

