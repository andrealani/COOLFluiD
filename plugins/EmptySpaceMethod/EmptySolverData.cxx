// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "EmptySpaceMethod/Empty.hh"
#include "EmptySpaceMethod/EmptySolverData.hh"
#include "EmptySpaceMethod/EmptyStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<
  NullMethodCommand< EmptySolverData >,EmptySolverData,EmptyModule >
  nullEmptySolverComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void EmptySolverData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("IntegratorOrder","Order of the Integration to be used for numerical quadrature.");
  options.addConfigOption< std::string >("IntegratorQuadrature","Type of Quadrature to be used in the Integration.");
  options.addConfigOption< std::string >("StrategyForSomething","A MethodStrategy to be used for some calculation (default = EmptyStrategy).");
}

//////////////////////////////////////////////////////////////////////////////

EmptySolverData::EmptySolverData(Common::SafePtr<Framework::Method> owner) :
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

  m_emptystrategyStr = "EmptyStrategy";
  setParameter( "StrategyForSomething", &m_emptystrategyStr );
}

//////////////////////////////////////////////////////////////////////////////

EmptySolverData::~EmptySolverData()
{
}

//////////////////////////////////////////////////////////////////////////////

void EmptySolverData::configure ( Config::ConfigArgs& args )
{
  SpaceMethodData::configure(args);
  SharedPtr< EmptySolverData > thisPtr(this);

  configureIntegrator();

  /* add here the setup for the specific integrators for each element that
   * have different set of shape and interpolator type */

  CFLog(INFO,"EmptySolver: integrator quadrature: " << m_intquadStr << "\n");
  CFLog(INFO,"EmptySolver: integrator order: " << m_intorderStr << "\n");

  /* add here different strategies configuration */

  CFLog(INFO,"Configure strategy type: " << m_emptystrategyStr << "\n");
  try {

    SafePtr< BaseMethodStrategyProvider< EmptySolverData, EmptyStrategy > >
      prov = Environment::Factory< EmptyStrategy >::getInstance().getProvider(
        m_emptystrategyStr );
    cf_assert(prov.isNotNull());
    m_emptystrategy = prov->create(m_emptystrategyStr,thisPtr);
    configureNested ( m_emptystrategy.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  EmptyStrategy ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< EmptySolverData, EmptyStrategy > >
      prov = Environment::Factory< EmptyStrategy >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_emptystrategy = prov->create("Null", thisPtr);

  }
  cf_assert(m_emptystrategy.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void EmptySolverData::configureIntegrator()
{
  const CFQuadrature::Type quadType = CFQuadrature::Convert::to_enum( m_intquadStr );
  const CFPolyOrder::Type order = CFPolyOrder::Convert::to_enum( m_intorderStr );

  m_volumeIntegrator.setIntegrationForAllGeo(quadType,order);
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<VolumeIntegrator> EmptySolverData::getVolumeIntegrator()
{
  return &m_volumeIntegrator;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

