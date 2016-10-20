// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/FluxReconstructionStrategy.hh"
#include "FluxReconstructionMethod/BaseInterfaceFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<
  NullMethodCommand< FluxReconstructionSolverData >,FluxReconstructionSolverData,FluxReconstructionModule >
  nullFluxReconstructionSolverComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolverData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("IntegratorOrder","Order of the Integration to be used for numerical quadrature.");
  options.addConfigOption< std::string >("IntegratorQuadrature","Type of Quadrature to be used in the Integration.");
  options.addConfigOption< std::string >("StrategyForSomething","A MethodStrategy to be used for some calculation (default = FluxReconstructionStrategy).");
  options.addConfigOption< std::string >("InterfaceFluxComputer","Name of the interface flux computer");
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolverData::FluxReconstructionSolverData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  m_lss(),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder()
{
  addConfigOptionsTo(this);

  m_intorderStr = "P1";
  m_intquadStr  = "INVALID";
  setParameter( "IntegratorOrder",      &m_intorderStr );
  setParameter( "IntegratorQuadrature", &m_intquadStr );

  m_fluxreconstructionstrategyStr = "FluxReconstructionStrategy";
  setParameter( "StrategyForSomething", &m_fluxreconstructionstrategyStr );
  
  m_interfacefluxStr = "BaseInterfaceFlux";
  setParameter( "InterfaceFluxComputer", &m_interfacefluxStr );
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolverData::~FluxReconstructionSolverData()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolverData::configure ( Config::ConfigArgs& args )
{
  SpaceMethodData::configure(args);
  SharedPtr< FluxReconstructionSolverData > thisPtr(this);

  configureIntegrator();

  /* add here the setup for the specific integrators for each element that
   * have different set of shape and interpolator type */

  CFLog(INFO,"FluxReconstructionSolver: integrator quadrature: " << m_intquadStr << "\n");
  CFLog(INFO,"FluxReconstructionSolver: integrator order: " << m_intorderStr << "\n");

  /* add here different strategies configuration */

  CFLog(INFO,"Configure strategy type: " << m_fluxreconstructionstrategyStr << "\n");
  try {

    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, FluxReconstructionStrategy > >
      prov = Environment::Factory< FluxReconstructionStrategy >::getInstance().getProvider(
        m_fluxreconstructionstrategyStr );
    cf_assert(prov.isNotNull());
    m_fluxreconstructionstrategy = prov->create(m_fluxreconstructionstrategyStr,thisPtr);
    configureNested ( m_fluxreconstructionstrategy.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  FluxReconstructionStrategy ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, FluxReconstructionStrategy > >
      prov = Environment::Factory< FluxReconstructionStrategy >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_fluxreconstructionstrategy = prov->create("Null", thisPtr);

  }
  cf_assert(m_fluxreconstructionstrategy.isNotNull());
  
  CFLog(INFO,"Configure strategy type: " << m_interfacefluxStr << "\n");
  try {

    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BaseInterfaceFlux > >
      prov = Environment::Factory< BaseInterfaceFlux >::getInstance().getProvider(
        m_interfacefluxStr );
    cf_assert(prov.isNotNull());
    m_interfaceflux = prov->create(m_interfacefluxStr,thisPtr);
    configureNested ( m_interfaceflux.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  BaseInterfaceFlux ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BaseInterfaceFlux > >
      prov = Environment::Factory< BaseInterfaceFlux >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_interfaceflux = prov->create("Null", thisPtr);

  }
  cf_assert(m_interfaceflux.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolverData::configureIntegrator()
{
  const CFQuadrature::Type quadType = CFQuadrature::Convert::to_enum( m_intquadStr );
  const CFPolyOrder::Type order = CFPolyOrder::Convert::to_enum( m_intorderStr );

  m_volumeIntegrator.setIntegrationForAllGeo(quadType,order);
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<VolumeIntegrator> FluxReconstructionSolverData::getVolumeIntegrator()
{
  return &m_volumeIntegrator;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

