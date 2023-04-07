// Copyright (C) 2016 KU Leuven, Belgium
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
#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/PrismFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/BasePointDistribution.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"
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
  options.addConfigOption< std::string >("LinearVar","Name of the linear variable set.");
  options.addConfigOption< std::string >("FluxPointDistribution","Name of the flux point distribution");
  options.addConfigOption< std::string >("RiemannFlux","Name of the Riemann flux.");
  options.addConfigOption< bool >("FreezeJacob","Flag telling whether to freeze the Jacobian.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("FreezeJacobIter","Iteration after which to freeze the Jacobian.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("FreezeJacobInterval","Amount of iterations to freeze the Jacobian before recalculation.");
  options.addConfigOption< CFreal >("DiffFluxDamping","Damping coefficient of diffusive flux scheme.");
  options.addConfigOption< bool >("AddArtificialViscosity","Flag telling whether to add artificial viscosity.");
  options.addConfigOption< std::string >("SolutionPointDistribution","Name of the solution point distribution");
  options.addConfigOption< std::string >("CorrectionFunctionComputer","Name of the correction function computer");
  options.addConfigOption< std::vector<std::string> >("BcTypes","Types of the boundary condition commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names of the boundary condition commands.");
  options.addConfigOption< bool >("ComputeVolumeForEachState" ,"Boolean telling whether to create a socket with the volume for each state, needed for some unsteady algorithms.");
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionSolverData::FluxReconstructionSolverData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  m_numJacob(CFNULL),
  m_lss(),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder(),
  m_faceBuilder(),
  m_faceBuilder2nd(),
  m_cellBuilder(),
  m_cellBuilder2nd(),
  m_frLocalData(),
  m_linearVarStr(),
  m_bcs(),
  m_bcsSP(),
  m_bcTypeStr(),
  m_bcNameStr(),
  m_bcTRSNameStr(),
  m_maxNbrRFluxPnts(),
  m_maxNbrStatesData(),
  m_resFactor(),
  m_hasDiffTerm(),
  m_bndFacesStartIdxs(),
  m_innerFacesStartIdxs(),
  m_partitionFacesStartIdxs(),
  m_updateToSolutionVecTrans(),
  m_solToUpdateInUpdateMatTrans()
{
  addConfigOptionsTo(this);
  
  m_linearVarStr = "Roe";
  setParameter( "LinearVar", &m_linearVarStr );
  
  m_riemannFluxStr = "RoeFlux";
  setParameter("RiemannFlux", &m_riemannFluxStr);
  
  m_diffDampCoeff = 1.0;
  setParameter("DiffFluxDamping", &m_diffDampCoeff);
  
  m_fluxPntDistributionStr = "Null";
  setParameter( "FluxPointDistribution", &m_fluxPntDistributionStr );
  
  m_solPntDistributionStr = "Null";
  setParameter( "SolutionPointDistribution", &m_solPntDistributionStr );
  
  m_correctionFunctionStr = "Null";
  setParameter( "CorrectionFunctionComputer", &m_correctionFunctionStr );
  
  m_freezeJacob = false;
  setParameter( "FreezeJacob", &m_freezeJacob );
  
  m_addAV = false;
  setParameter( "AddArtificialViscosity", &m_addAV );

  m_createVolumesSocketBool = true;
  setParameter("ComputeVolumeForEachState", &m_createVolumesSocketBool);
  
  // options for bc commands
  m_bcTypeStr = vector<std::string>();
  setParameter("BcTypes",&m_bcTypeStr);

  m_bcNameStr = vector<std::string>();
  setParameter("BcNames",&m_bcNameStr);
  
  m_freezeJacobIter = MathTools::MathConsts::CFuintMax();
  setParameter( "FreezeJacobIter", &m_freezeJacobIter );
  
  m_freezeJacobInterval = 1;
  setParameter( "FreezeJacobInterval", &m_freezeJacobInterval );
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
  m_faceBuilder2nd.setup();
  
  // setup cell to face builder
  m_cellBuilder.setup();
  m_cellBuilder2nd.setup();
  
  // create numerical Jacobian computer
  m_numJacob.reset(new NumericalJacobian("NumericalJacobian"));
  
  // set reference values in numerical Jacobian computer
  RealVector refValues = PhysicalModelStack::getActive()->getImplementor()->getRefStateValues();
  m_numJacob->setRefValues(refValues);
  
  // setup the variable sets
  _updateVar   ->setup();
  _solutionVar ->setup();
  _diffusiveVar->setup();
  
  // create local FR data
  createFRLocalData();
  
  // setup StatesReconstructor
  m_statesReconstructor->setup();
  
  // set the hasDiffTerm boolean
  /// @note it would be better to check a name related to the DiffusiveVarSet here
  m_hasDiffTerm = (_diffusiveVarStr != "Null");
  
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
  
  std::string solToUpdateInUpdateMatTransStr =
    VarSetMatrixTransformer::getProviderName
    (physModel->getConvectiveName(), _solutionVarStr,
     _updateVarStr, _updateVarStr);

  CFLog(VERBOSE, "Configuring VarSet Transformer: " <<
	solToUpdateInUpdateMatTransStr << "\n");
  
  Common::SafePtr<VarSetMatrixTransformer::PROVIDER> matTransProv = CFNULL;

  try {
    matTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetMatrixTransformer, solToUpdateInUpdateMatTransStr);
  }
  catch (Common::NoSuchValueException& e) {
    solToUpdateInUpdateMatTransStr = "Identity";

    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing IdentityVarSetMatrixTransformer instead ..." << "\n");
    matTransProv = FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetMatrixTransformer, solToUpdateInUpdateMatTransStr);
  }

  cf_assert(matTransProv.isNotNull());
  m_solToUpdateInUpdateMatTrans.reset(matTransProv->create(physModel->getImplementor()));
  cf_assert(m_solToUpdateInUpdateMatTrans.getPtr() != CFNULL);
  
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
  
  // Configure flux point distribution
  CFLog(INFO,"Configure flux point distribution: " << m_fluxPntDistributionStr << "\n");
  try {

    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BasePointDistribution > >
      prov = Environment::Factory< BasePointDistribution >::getInstance().getProvider(
        m_fluxPntDistributionStr );
    cf_assert(prov.isNotNull());
    m_fluxPntDistribution = prov->create(m_fluxPntDistributionStr,thisPtr);
    configureNested ( m_fluxPntDistribution.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  BasePointDistribution ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BasePointDistribution > >
      prov = Environment::Factory< BasePointDistribution >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_fluxPntDistribution = prov->create("Null", thisPtr);

  }
  cf_assert(m_fluxPntDistribution.isNotNull());
  
  // Configure solution point distribution
  CFLog(INFO,"Configure solution point distribution: " << m_solPntDistributionStr << "\n");
  try {

    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BasePointDistribution > >
      prov = Environment::Factory< BasePointDistribution >::getInstance().getProvider(
        m_solPntDistributionStr );
    cf_assert(prov.isNotNull());
    m_solPntDistribution = prov->create(m_solPntDistributionStr,thisPtr);
    configureNested ( m_solPntDistribution.getPtr(), args );

  } catch (Common::NoSuchValueException& e) {

    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  BasePointDistribution ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BasePointDistribution > >
      prov = Environment::Factory< BasePointDistribution >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_solPntDistribution = prov->create("Null", thisPtr);

  }
  cf_assert(m_solPntDistribution.isNotNull());

  CFLog(INFO,"Configure strategy type: " << m_correctionFunctionStr << "\n");
  try {
        
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BaseCorrectionFunction > >
      prov = Environment::Factory< BaseCorrectionFunction >::getInstance().getProvider(
        m_correctionFunctionStr );
    cf_assert(prov.isNotNull());
    m_correctionFunction = prov->create(m_correctionFunctionStr,thisPtr);
    configureNested ( m_correctionFunction.getPtr(), args );
        
  } catch (Common::NoSuchValueException& e) {
        
    CFLog(INFO, e.what() << "\n");
    CFLog(INFO, "Choosing Null of type: " <<  BaseCorrectionFunction ::getClassName() << " instead...\n");
    SafePtr< BaseMethodStrategyProvider< FluxReconstructionSolverData, BaseCorrectionFunction > >
      prov = Environment::Factory< BaseCorrectionFunction >::getInstance().getProvider("Null");
    cf_assert(prov.isNotNull());
    m_correctionFunction = prov->create("Null", thisPtr);
        
  }
  cf_assert(m_correctionFunction.isNotNull());
  
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
  
  // Configure states reconstructor
  Common::SafePtr<BaseMethodStrategyProvider< FluxReconstructionSolverData , ReconstructStatesFluxReconstruction > > recProv =
    Environment::Factory< ReconstructStatesFluxReconstruction >
                                                  ::getInstance().getProvider("ReconstructStatesFluxReconstruction");
  cf_assert(recProv.isNotNull());
  m_statesReconstructor = recProv->create("ReconstructStatesFluxReconstruction",thisPtr);

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
        m_frLocalData[iElemType] = new TriagFluxReconstructionElementData(polyOrder,getSolPntDistribution(),getFluxPntDistribution());
      } break;
      case CFGeoShape::QUAD:
      {
        m_frLocalData[iElemType] = new QuadFluxReconstructionElementData(polyOrder,getSolPntDistribution(),getFluxPntDistribution());
      } break;
      case CFGeoShape::TETRA:
      {
        m_frLocalData[iElemType] = new TetraFluxReconstructionElementData(polyOrder,getSolPntDistribution(),getFluxPntDistribution());      
      } break;
      case CFGeoShape::HEXA:
      {
        m_frLocalData[iElemType] = new HexaFluxReconstructionElementData(polyOrder,getSolPntDistribution(),getFluxPntDistribution());
      } break;
      case CFGeoShape::PRISM:
      {
        m_frLocalData[iElemType] = new PrismFluxReconstructionElementData(polyOrder,getSolPntDistribution(),getFluxPntDistribution());      
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

