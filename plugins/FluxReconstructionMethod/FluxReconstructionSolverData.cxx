// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/FilesystemException.hh"
#include "Framework/NamespaceSwitcher.hh"

//#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/QuadFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/HexaFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/BasePointDistribution.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/BaseBndFaceTermComputer.hh"
#include "FluxReconstructionMethod/BaseFaceTermComputer.hh"
#include "FluxReconstructionMethod/BaseVolTermComputer.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/ReconstructStatesFluxReconstruction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

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
  options.addConfigOption< std::string >("LinearVar","Name of the linear variable set.");
  options.addConfigOption< std::string >("FluxPointDistribution","Name of the flux point distribution");
  options.addConfigOption< std::string >("SolutionPointDistribution","Name of the solution point distribution");
  options.addConfigOption< std::string >("RiemannFlux","Name of the Riemann flux.");
  options.addConfigOption< std::string >("BndFaceTermComputer","Name of the boundary face term computer."  );
  options.addConfigOption< std::string >("FaceTermComputer","Name of the face term computer."  );
  options.addConfigOption< std::string >("VolTermComputer" ,"Name of the volume term computer.");
  options.addConfigOption< std::string >("CorrectionFunctionComputer","Name of the correction function computer");
  options.addConfigOption< std::vector<std::string> >("BcTypes","Types of the boundary condition commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names of the boundary condition commands.");
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolverData::FluxReconstructionSolverData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  m_lss(),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder(),
  m_faceBuilder(),
  m_statesReconstructor(),
  m_volTermComputerStr(),
  m_volTermComputer(),
  m_faceTermComputerStr(),
  m_faceTermComputer(),
  m_bndFaceTermComputerStr(),
  m_bndFaceTermComputer(),
  m_linearVarStr(),
  m_riemannFluxStr(),
  m_riemannFlux(),
  m_bcs(),
  m_bcsSP(),
  m_bcTypeStr(),
  m_bcNameStr(),
  m_bcTRSNameStr(),
  m_innerFacesStartIdxs(),
  m_bndFacesStartIdxs(),
  m_frLocalData(),
  m_maxNbrStatesData(),
  m_maxNbrRFluxPnts(),
  m_resFactor(),
  m_updateToSolutionVecTrans()
  //socket_solCoords1D("solCoords1D"),
  //socket_flxCoords1D("flxCoords1D")
{
  addConfigOptionsTo(this);

  //m_intorderStr = "P1";
  //m_intquadStr  = "INVALID";
  //setParameter( "IntegratorOrder",      &m_intorderStr );
  //setParameter( "IntegratorQuadrature", &m_intquadStr );
  
  m_linearVarStr = "Roe";
  setParameter( "LinearVar", &m_linearVarStr );
  
  m_riemannFluxStr = "RoeFlux";
  setParameter("RiemannFlux", &m_riemannFluxStr);
  
  m_fluxpntdistributionStr = "Null";
  setParameter( "FluxPointDistribution", &m_fluxpntdistributionStr );
  
  m_solpntdistributionStr = "Null";
  setParameter( "SolutionPointDistribution", &m_solpntdistributionStr );
  
  m_correctionfunctionStr = "Null";
  setParameter( "CorrectionFunctionComputer", &m_correctionfunctionStr );
  
  m_bndFaceTermComputerStr = "BaseBndFaceTermComputer";
  setParameter("BndFaceTermComputer", &m_bndFaceTermComputerStr);

  m_faceTermComputerStr = "BaseFaceTermComputer";
  setParameter("FaceTermComputer", &m_faceTermComputerStr);

  m_volTermComputerStr = "BaseVolTermComputer";
  setParameter("VolTermComputer", &m_volTermComputerStr);
  
  // options for bc commands
  m_bcTypeStr = vector<std::string>();
  setParameter("BcTypes",&m_bcTypeStr);

  m_bcNameStr = vector<std::string>();
  setParameter("BcNames",&m_bcNameStr);
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolverData::~FluxReconstructionSolverData()
{
}

//////////////////////////////////////////////////////////////////////////////

// std::vector< Common::SafePtr< BaseDataSocketSource > >
//   FluxReconstructionSolverData::providesSockets()
// {
//   std::vector< Common::SafePtr< BaseDataSocketSource > > result;
//   result.push_back(&socket_solCoords1D);
//   result.push_back(&socket_flxCoords1D);
//   return result;
// }

//////////////////////////////////////////////////////////////////////////////



void FluxReconstructionSolverData::setup()
{
  CFAUTOTRACE;
  
  SpaceMethodData::setup();
  
  // setup TRS Geo builder
  m_stdTrsGeoBuilder.setup();
  
  // setup face builder
  m_faceBuilder.setup();
  
  // setup the variable sets
  _updateVar   ->setup();
  _solutionVar ->setup();
  _diffusiveVar->setup();
  
  // create local FR data
  createFRLocalData();
  
  // setup StatesReconstructor
  m_statesReconstructor->setup();
  
  // create the transformer from update to solution variables 
  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  Common::SafePtr<VarSetTransformer::PROVIDER> vecTransProv = CFNULL;
  
  std::string updateToSolutionVecTransStr =
  VarSetTransformer::getProviderName
  (physModel->getConvectiveName(), _updateVarStr, _solutionVarStr);
  
  CFLog(VERBOSE, "Configuring VarSet Transformer: " <<
        updateToSolutionVecTransStr << "\n");
  
  try {
    vecTransProv = Environment::Factory<VarSetTransformer>::getInstance().getProvider
    (updateToSolutionVecTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    updateToSolutionVecTransStr = "Identity";
    
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetTransformer instead ..." << "\n");
    vecTransProv = Environment::Factory<VarSetTransformer>::getInstance().getProvider
    (updateToSolutionVecTransStr);
  }
  CFLog(VERBOSE, "FluxReconstructionSolverData::setup() => updateToSolutionVarName = " << updateToSolutionVecTransStr << "\n");
  
  cf_assert(vecTransProv.isNotNull());
  m_updateToSolutionVecTrans.reset(vecTransProv->create(physModel->getImplementor()));
  cf_assert(m_updateToSolutionVecTrans.isNotNull());
  m_updateToSolutionVecTrans->setup(2);
  
//   // setup socket solCoords1D
//   DataHandle<std::vector<CFreal> > solCoords1D = socket_solCoords1D.getDataHandle();
//   SafePtr< std::vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
//   // get the order of the polynomial interpolation
//   const CFPolyOrder::Type polyOrder = static_cast<CFPolyOrder::Type>((*elemType)[0].getSolOrder());
//   solCoords1D.resize(polyOrder+1);
//   solCoords1D = getSolPntDistribution()->getLocalCoords1D(polyOrder);
//   
//   // setup socket flxCoords1D
//   DataHandle<std::vector<CFreal> > flxCoords1D = socket_flxCoords1D.getDataHandle();
//   flxCoords1D.resize(polyOrder+1);
//   flxCoords1D = getFluxPntDistribution()->getLocalCoords1D(polyOrder);
  
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
  
  // Configure boundary face term computer
  CFLog(INFO,"Configure boundary face term computer: " << m_bndFaceTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , BaseBndFaceTermComputer > > prov =
      Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider(m_bndFaceTermComputerStr);
    cf_assert(prov.isNotNull());
    m_bndFaceTermComputer = prov->create(m_bndFaceTermComputerStr,thisPtr);
    configureNested ( m_bndFaceTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing BaseBndFaceTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , BaseBndFaceTermComputer > > prov =
      Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider("BaseBndFaceTermComputer");
    cf_assert(prov.isNotNull());
    m_bndFaceTermComputer = prov->create("BaseBndFaceTermComputer", thisPtr);
    configureNested ( m_bndFaceTermComputer.getPtr(), args );
  }
  cf_assert(m_bndFaceTermComputer.isNotNull());
  
  // Configure face term computer
  CFLog(INFO,"Configure face term computer: " << m_faceTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , BaseFaceTermComputer > > prov =
      Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider(m_faceTermComputerStr);
    cf_assert(prov.isNotNull());
    m_faceTermComputer = prov->create(m_faceTermComputerStr,thisPtr);
    configureNested ( m_faceTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing BaseFaceTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , BaseFaceTermComputer > > prov =
      Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider("BaseFaceTermComputer");
    cf_assert(prov.isNotNull());
    m_faceTermComputer = prov->create("BaseFaceTermComputer", thisPtr);
    configureNested ( m_faceTermComputer.getPtr(), args );
  }
  cf_assert(m_faceTermComputer.isNotNull());
  
  // Configure volume term computer
  CFLog(INFO,"Configure volume term computer: " << m_volTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider(m_volTermComputerStr);
    cf_assert(prov.isNotNull());
    m_volTermComputer = prov->create(m_volTermComputerStr,thisPtr);
    configureNested ( m_volTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing BaseVolTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider("BaseVolTermComputer");
    cf_assert(prov.isNotNull());
    m_volTermComputer = prov->create("BaseVolTermComputer", thisPtr);
    configureNested ( m_volTermComputer.getPtr(), args );
  }
  cf_assert(m_volTermComputer.isNotNull());
  
  // Configure flux point distribution
  CFLog(INFO,"Configure flux point distribution: " << m_fluxpntdistributionStr << "\n");
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
  
  // Configure solution point distribution
  CFLog(INFO,"Configure solution point distribution: " << m_solpntdistributionStr << "\n");
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

  CFLog(INFO,"Configure strategy type: " << m_correctionfunctionStr << "\n");
  try {
        
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BaseCorrectionFunction > >
      prov = Environment::Factory< BaseCorrectionFunction >::getInstance().getProvider(
        m_correctionfunctionStr );
    cf_assert(prov.isNotNull());
    m_correctionfunction = prov->create(m_correctionfunctionStr,thisPtr);
    configureNested ( m_correctionfunction.getPtr(), args );
        
  } catch (Common::NoSuchValueException& e) {
        
    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  BaseCorrectionFunction ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BaseCorrectionFunction > >
      prov = Environment::Factory< BaseCorrectionFunction >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_correctionfunction = prov->create("Null", thisPtr);
        
  }
  cf_assert(m_correctionfunction.isNotNull());
  
  // Configure Riemann flux
  CFLog(INFO,"Configure Riemann flux: " << m_riemannFluxStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , RiemannFlux > > prov =
      Environment::Factory< RiemannFlux >::getInstance().getProvider(m_riemannFluxStr);
    cf_assert(prov.isNotNull());
    m_riemannFlux = prov->create(m_riemannFluxStr,thisPtr);
    configureNested ( m_riemannFlux.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing RoeFlux instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , RiemannFlux > > prov =
      Environment::Factory< RiemannFlux >::getInstance().getProvider("RoeFlux");
    cf_assert(prov.isNotNull());
    m_riemannFlux = prov->create("RoeFlux", thisPtr);
    configureNested ( m_riemannFlux.getPtr(), args );
  }
  cf_assert(m_riemannFlux.isNotNull());
  
  // Configure BC state computers
  CFLog(INFO,"Configure BC state computers\n");
  cf_assert(m_bcTypeStr.size() == m_bcNameStr.size());

  // number of boundary conditions
  const CFuint nbrBcs = m_bcTypeStr.size();

  // resize m_bcs and m_bcsSP
  m_bcs.resize(nbrBcs);
  m_bcsSP.resize(nbrBcs);

  for (CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    CFLog(INFO, "BC type = " << m_bcTypeStr[iBc] << "\n");
    CFLog(INFO, "BC name = " << m_bcNameStr[iBc] << "\n");

    Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , BCStateComputer > > prov =
      Environment::Factory<BCStateComputer>::getInstance().getProvider(m_bcTypeStr[iBc]);
    cf_assert(prov.isNotNull());
    m_bcs[iBc] = prov->create(m_bcNameStr[iBc],thisPtr);
    configureNested(m_bcs[iBc].getPtr(), args);
    cf_assert(m_bcs[iBc].isNotNull());

    // set SafePtr
    m_bcsSP[iBc] = m_bcs[iBc].getPtr();
  }
  
  // Configure states reconstructor
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

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

