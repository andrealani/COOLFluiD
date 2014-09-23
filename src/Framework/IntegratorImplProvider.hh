// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IntegratorImplProvider_hh
#define COOLFluiD_Framework_IntegratorImplProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"
#include "Common/SelfRegistPtr.hh"

#include "Environment/ConcreteProvider.hh"
#include "Framework/IntegratorRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a factory for IntegratorImpl's.
/// @author Tiago Quintino
/// @author Andrea Lani
template <typename BASE_INTEGRATOR, typename INTEGRATOR>
class IntegratorImplProvider : public Environment::ConcreteProvider<BASE_INTEGRATOR> {
public:

  /// Constructor.
  /// @see ConcreteProvider
  IntegratorImplProvider() :
    Environment::ConcreteProvider<BASE_INTEGRATOR>(INTEGRATOR::getName())
  {
    IntegratorRegister<BASE_INTEGRATOR>::getInstance().
    registIntegrator(
      INTEGRATOR::getName(),
      INTEGRATOR::getIntegrationType(),
      INTEGRATOR::getQuadratureType(),
      INTEGRATOR::getIntegrationOrder(),
      INTEGRATOR::getShape(),
      INTEGRATOR::getInterpolatorType(),
      INTEGRATOR::getInterpolatorOrder());
  }

  /// Default destructor
  ~IntegratorImplProvider() {}

  /// Create the required object of dynamical type = BASE_INTEGRATOR
  Common::SelfRegistPtr<BASE_INTEGRATOR> create()
  {
    return Common::SelfRegistPtr<BASE_INTEGRATOR>(new INTEGRATOR(), this);
  }

  /// Free an instance created by this factory
  /// @param ptr pointer to be freed
  void freeInstance (void * ptr)
  {
    cf_assert(ptr != CFNULL);

    INTEGRATOR* obj = reinterpret_cast<INTEGRATOR*>(ptr);
    cf_assert(obj != CFNULL);
    deletePtr<INTEGRATOR>( obj );
  }

}; // end of class IntegratorImplProvider

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IntegratorImplProvider_hh
