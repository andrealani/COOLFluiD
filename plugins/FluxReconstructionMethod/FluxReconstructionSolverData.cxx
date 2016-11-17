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
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/QuadFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionStrategy.hh"
#include "FluxReconstructionMethod/BaseInterfaceFlux.hh"
#include "FluxReconstructionMethod/BaseFluxPntDistribution.hh"

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
  options.addConfigOption< std::string >("FluxPointDistribution","Name of the flux point distribution");
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolverData::FluxReconstructionSolverData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  m_lss(),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder(),
  m_frLocalData()
{
  addConfigOptionsTo(this);

  m_intorderStr = "P1";
  m_intquadStr  = "INVALID";
  setParameter( "IntegratorOrder",      &m_intorderStr );
  setParameter( "IntegratorQuadrature", &m_intquadStr );

  m_fluxreconstructionstrategyStr = "FluxReconstructionStrategy";
  setParameter( "StrategyForSomething", &m_fluxreconstructionstrategyStr );
  
  m_interfacefluxStr = "Null";
  setParameter( "InterfaceFluxComputer", &m_interfacefluxStr );
  
  m_fluxpntdistributionStr = "Null";
  setParameter( "FluxPointDistribution", &m_fluxpntdistributionStr );
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolverData::~FluxReconstructionSolverData()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolverData::setup()
{
  CFAUTOTRACE;
  
  SpaceMethodData::setup();
  
  // create local FR data
  createFRLocalData();
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
  
  CFLog(INFO,"Configure strategy type: " << m_fluxpntdistributionStr << "\n");
  try {

    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BaseFluxPntDistribution > >
      prov = Environment::Factory< BaseFluxPntDistribution >::getInstance().getProvider(
        m_fluxpntdistributionStr );
    cf_assert(prov.isNotNull());
    m_fluxpntdistribution = prov->create(m_fluxpntdistributionStr,thisPtr);
    configureNested ( m_fluxpntdistribution.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  BaseFluxPntDistribution ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BaseFluxPntDistribution > >
      prov = Environment::Factory< BaseFluxPntDistribution >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_fluxpntdistribution = prov->create("Null", thisPtr);

  }
  cf_assert(m_fluxpntdistribution.isNotNull());

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

void FluxReconstructionSolverData::createFRLocalData()
{
  CFAUTOTRACE;

  // get the ElementTypeData
  SafePtr< std::vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // resize frLocalData
  m_frLocalData.resize(nbrElemTypes);

  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get the element shape
    const CFGeoShape::Type elemShape = (*elemType)[iElemType].getGeoShape();

    // get the order of the polynomial interpolation
    const CFPolyOrder::Type polyOrder = static_cast<CFPolyOrder::Type>((*elemType)[iElemType].getSolOrder());

    switch(elemShape)
    {
      case CFGeoShape::LINE:
      {
        throw Common::NotImplementedException (FromHere(),"FR has not been implemented for 1D");
      } break;
      case CFGeoShape::TRIAG:
      {
        throw Common::NotImplementedException (FromHere(),"FR has not been implemented for triangular cells");
      } break;
      case CFGeoShape::QUAD:
      {
        m_frLocalData[iElemType] = new QuadFluxReconstructionElementData(polyOrder);
      } break;
      case CFGeoShape::TETRA:
      {
        throw Common::NotImplementedException (FromHere(),"FR has not been implemented for tetrahedral cells");
      } break;
      //case CFGeoShape::HEXA:
      //{
      //  m_frLocalData[iElemType] = new HexaSpectralFDElementData(polyOrder);
      //} break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"FR method not implemented for elements of type "
                                      + StringOps::to_str(elemShape) + ".");
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< FluxReconstructionElementData* >& FluxReconstructionSolverData::getFRLocalData()
{
  return m_frLocalData;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

