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
#include "FluxReconstructionMethod/HexaFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/BaseInterfaceFlux.hh"
#include "FluxReconstructionMethod/BasePointDistribution.hh"
#include "FluxReconstructionMethod/ReconstructStatesFluxReconstruction.hh"

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
  //options.addConfigOption< std::string >("IntegratorOrder","Order of the Integration to be used for numerical quadrature.");
  //options.addConfigOption< std::string >("IntegratorQuadrature","Type of Quadrature to be used in the Integration.");
  options.addConfigOption< std::string >("InterfaceFluxComputer","Name of the interface flux computer");
  options.addConfigOption< std::string >("FluxPointDistribution","Name of the flux point distribution");
  options.addConfigOption< std::string >("SolutionPointDistribution","Name of the solution point distribution");
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolverData::FluxReconstructionSolverData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  m_lss(),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder(),
  m_frLocalData(),
  m_statesReconstructor(),
  m_faceBuilder()
{
  addConfigOptionsTo(this);

  //m_intorderStr = "P1";
  //m_intquadStr  = "INVALID";
  //setParameter( "IntegratorOrder",      &m_intorderStr );
  //setParameter( "IntegratorQuadrature", &m_intquadStr );
  
  m_interfacefluxStr = "Null";
  setParameter( "InterfaceFluxComputer", &m_interfacefluxStr );
  
  m_fluxpntdistributionStr = "Null";
  setParameter( "FluxPointDistribution", &m_fluxpntdistributionStr );
  
  m_solpntdistributionStr = "Null";
  setParameter( "SolutionPointDistribution", &m_solpntdistributionStr );
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
  
  // setup TRS Geo builder
  m_stdTrsGeoBuilder.setup();
  
  // setup face builder
  m_faceBuilder.setup();
  
  // create local FR data
  createFRLocalData();
  
  // setup StatesReconstructor
  m_statesReconstructor->setup();
  
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolverData::unsetup()
{
  CFAUTOTRACE;
  
  SpaceMethodData::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionSolverData::configure ( Config::ConfigArgs& args )
{
  SpaceMethodData::configure(args);
  SharedPtr< FluxReconstructionSolverData > thisPtr(this);

  //configureIntegrator();

//   /* add here the setup for the specific integrators for each element that
//    * have different set of shape and interpolator type */
// 
//   CFLog(INFO,"FluxReconstructionSolver: integrator quadrature: " << m_intquadStr << "\n");
//   CFLog(INFO,"FluxReconstructionSolver: integrator order: " << m_intorderStr << "\n");

  /* add here different strategies configuration */
  
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

    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BasePointDistribution > >
      prov = Environment::Factory< BasePointDistribution >::getInstance().getProvider(
        m_fluxpntdistributionStr );
    cf_assert(prov.isNotNull());
    m_fluxpntdistribution = prov->create(m_fluxpntdistributionStr,thisPtr);
    configureNested ( m_fluxpntdistribution.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  BasePointDistribution ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BasePointDistribution > >
      prov = Environment::Factory< BasePointDistribution >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_fluxpntdistribution = prov->create("Null", thisPtr);

  }
  cf_assert(m_fluxpntdistribution.isNotNull());
  
  CFLog(INFO,"Configure strategy type: " << m_solpntdistributionStr << "\n");
  try {

    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BasePointDistribution > >
      prov = Environment::Factory< BasePointDistribution >::getInstance().getProvider(
        m_solpntdistributionStr );
    cf_assert(prov.isNotNull());
    m_solpntdistribution = prov->create(m_solpntdistributionStr,thisPtr);
    configureNested ( m_solpntdistribution.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  BasePointDistribution ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BasePointDistribution > >
      prov = Environment::Factory< BasePointDistribution >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_solpntdistribution = prov->create("Null", thisPtr);

  }
  cf_assert(m_solpntdistribution.isNotNull());
  
  // states reconstructor
  Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , ReconstructStatesFluxReconstruction > > recProv =
    Environment::Factory< ReconstructStatesFluxReconstruction >
                                                  ::getInstance().getProvider("ReconstructStatesFluxReconstruction");
  cf_assert(recProv.isNotNull());
  m_statesReconstructor = recProv->create("ReconstructStatesFluxReconstruction",thisPtr);

}

//////////////////////////////////////////////////////////////////////////////

// void FluxReconstructionSolverData::configureIntegrator()
// {
//   const CFQuadrature::Type quadType = CFQuadrature::Convert::to_enum( m_intquadStr );
//   const CFPolyOrder::Type order = CFPolyOrder::Convert::to_enum( m_intorderStr );
// 
//   m_volumeIntegrator.setIntegrationForAllGeo(quadType,order);
// }

//////////////////////////////////////////////////////////////////////////////

// SafePtr<VolumeIntegrator> FluxReconstructionSolverData::getVolumeIntegrator()
// {
//   return &m_volumeIntegrator;
// }

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
        m_frLocalData[iElemType] = new TriagFluxReconstructionElementData(polyOrder,getSolPntDistribution(),getFluxPntDistribution());
      } break;
      case CFGeoShape::QUAD:
      {
        m_frLocalData[iElemType] = new QuadFluxReconstructionElementData(polyOrder,getSolPntDistribution(),getFluxPntDistribution());
      } break;
      case CFGeoShape::TETRA:
      {
        throw Common::NotImplementedException (FromHere(),"FR has not been implemented for tetrahedral cells");
      } break;
      case CFGeoShape::HEXA:
      {
        m_frLocalData[iElemType] = new HexaFluxReconstructionElementData(polyOrder,getSolPntDistribution(),getFluxPntDistribution());
      } break;
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

SafePtr< ReconstructStatesFluxReconstruction > FluxReconstructionSolverData::getStatesReconstructor()
// this function has to be put in the implementation file because of the forward declaration of ReconstructStatesSpectralFD
{
  return m_statesReconstructor.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

