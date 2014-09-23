#include "Common/FilesystemException.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/BaseFaceTermComputer.hh"
#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/FaceDiffusiveFlux.hh"
#include "SpectralFD/HexaSpectralFDElementData.hh"
#include "SpectralFD/QuadSpectralFDElementData.hh"
#include "SpectralFD/ReconstructStatesSpectralFD.hh"
#include "SpectralFD/RiemannFlux.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<
  NullMethodCommand< SpectralFDMethodData >,SpectralFDMethodData,SpectralFDModule >
  nullSpectralFDMethodComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("LinearVar"       ,"Name of the linear variable set." );
  options.addConfigOption< std::string >("RiemannFlux"     ,"Name of the Riemann flux."        );
  options.addConfigOption< std::string >("FaceDiffFlux"    ,"Name of the face diffusive flux." );
  options.addConfigOption< std::string >("BndFaceTermComputer","Name of the boundary face term computer."  );
  options.addConfigOption< std::string >("FaceTermComputer","Name of the face term computer."  );
  options.addConfigOption< std::string >("VolTermComputer" ,"Name of the volume term computer.");
  options.addConfigOption< std::vector<std::string> >("BcTypes","Types of the boundary condition commands.");
  options.addConfigOption< std::vector<std::string> >("BcNames","Names of the boundary condition commands.");
  options.addConfigOption< bool >("ComputeVolumeForEachState" ,"Boolean telling whether to create a socket with the volume for each state, needed for some unsteady algorithms.");
  options.addConfigOption< std::string >("InterpolationType","string defining the interpolation type to use (standard or optimized)");
}

//////////////////////////////////////////////////////////////////////////////

SpectralFDMethodData::SpectralFDMethodData(Common::SafePtr<Framework::Method> owner) :
  SpaceMethodData(owner),
  m_numJacob(CFNULL),
  m_convergenceMtd(),
  m_stdTrsGeoBuilder(),
  m_faceBuilder(),
  m_cellBuilder(),
  m_cellBuilder2nd(),
  m_statesReconstructor(),
  m_volTermComputerStr(),
  m_volTermComputer(),
  m_volTermComputer2nd(),
  m_faceTermComputerStr(),
  m_faceTermComputer(),
  m_addFaceTermComputers(),
  m_addFaceTermComputersSP(),
  m_bndFaceTermComputerStr(),
  m_bndFaceTermComputer(),
  m_addBndFaceTermComputers(),
  m_addBndFaceTermComputersSP(),
  m_linearVarStr(),
  m_riemannFluxStr(),
  m_faceDiffFluxStr(),
  m_riemannFlux(),
  m_faceDiffFlux(),
  m_bcs(),
  m_bcsSP(),
  m_bcTypeStr(),
  m_bcNameStr(),
  m_bcTRSNameStr(),
  m_innerFacesStartIdxs(),
  m_bndFacesStartIdxs(),
  m_sdLocalData(),
  m_maxNbrStatesData(),
  m_maxNbrRFluxPnts(),
  m_hasDiffTerm(),
  m_resFactor(),
  m_createVolumesSocketBool(),
  m_interpolationType(),
  m_3StepsTMSparams(),
  m_updateToSolutionVecTrans()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_linearVarStr = "Roe";
  setParameter( "LinearVar", &m_linearVarStr );

  m_riemannFluxStr = "RoeFlux";
  setParameter("RiemannFlux", &m_riemannFluxStr);

  m_faceDiffFluxStr = "LocalApproach";
  setParameter("FaceDiffFlux", &m_faceDiffFluxStr);

  m_bndFaceTermComputerStr = "BaseBndFaceTermComputer";
  setParameter("BndFaceTermComputer", &m_bndFaceTermComputerStr);

  m_faceTermComputerStr = "BaseFaceTermComputer";
  setParameter("FaceTermComputer", &m_faceTermComputerStr);

  m_volTermComputerStr = "BaseVolTermComputer";
  setParameter("VolTermComputer", &m_volTermComputerStr);

  m_createVolumesSocketBool = false;
  setParameter("ComputeVolumeForEachState", &m_createVolumesSocketBool);

  m_interpolationType = "standard";
  setParameter("InterpolationType", &m_interpolationType);

  // options for bc commands
  m_bcTypeStr = vector<std::string>();
  setParameter("BcTypes",&m_bcTypeStr);

  m_bcNameStr = vector<std::string>();
  setParameter("BcNames",&m_bcNameStr);
}

//////////////////////////////////////////////////////////////////////////////

SpectralFDMethodData::~SpectralFDMethodData()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::configure ( Config::ConfigArgs& args )
{
  SpaceMethodData::configure(args);
  SharedPtr< SpectralFDMethodData > thisPtr(this);

  /* add here different strategies configuration */
  CFLog(INFO,"SpectralFDMethod: configureBndFaceTermComputer()\n");
  configureBndFaceTermComputer(args);

  CFLog(INFO,"SpectralFDMethod: configureFaceTermComputer()\n");
  configureFaceTermComputer(args);

  CFLog(INFO,"SpectralFDMethod: configureVolTermComputer()\n");
  configureVolTermComputer(args);

  CFLog(INFO,"SpectralFDMethod: configureBCStateComputers()\n");
  configureBCStateComputers(args);

  CFLog(INFO,"SpectralFDMethod: configureRiemannFlux()\n");
  configureRiemannFlux(args);

  CFLog(INFO,"SpectralFDMethod: configureFaceDiffusiveFlux()\n");
  configureFaceDiffusiveFlux(args);

  // states reconstructor
  Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , ReconstructStatesSpectralFD > > recProv =
    Environment::Factory< ReconstructStatesSpectralFD >
                                                  ::getInstance().getProvider("ReconstructStatesSpectralFD");
  cf_assert(recProv.isNotNull());
  m_statesReconstructor = recProv->create("ReconstructStatesSpectralFD",thisPtr);
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::configureRiemannFlux( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFDMethodData> thisPtr(this);

  CFLogDebugMin("SpectralFD: Using Riemann flux: " << m_riemannFluxStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , RiemannFlux > > prov =
      Environment::Factory< RiemannFlux >::getInstance().getProvider(m_riemannFluxStr);
    cf_assert(prov.isNotNull());
    m_riemannFlux = prov->create(m_riemannFluxStr,thisPtr);
    configureNested ( m_riemannFlux.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing RoeFlux instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , RiemannFlux > > prov =
      Environment::Factory< RiemannFlux >::getInstance().getProvider("RoeFlux");
    cf_assert(prov.isNotNull());
    m_riemannFlux = prov->create("RoeFlux", thisPtr);
    configureNested ( m_riemannFlux.getPtr(), args );
  }
  cf_assert(m_riemannFlux.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::configureFaceDiffusiveFlux( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFDMethodData> thisPtr(this);

  CFLogDebugMin("SpectralFD: Using face diffusive flux: " << m_faceDiffFluxStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , FaceDiffusiveFlux > > prov =
      Environment::Factory< FaceDiffusiveFlux >::getInstance().getProvider(m_faceDiffFluxStr);
    cf_assert(prov.isNotNull());
    m_faceDiffFlux = prov->create(m_faceDiffFluxStr,thisPtr);
    configureNested ( m_faceDiffFlux.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing RoeFlux instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , FaceDiffusiveFlux > > prov =
      Environment::Factory< FaceDiffusiveFlux >::getInstance().getProvider("LocalApproach");
    cf_assert(prov.isNotNull());
    m_faceDiffFlux = prov->create("LocalApproach", thisPtr);
    configureNested ( m_faceDiffFlux.getPtr(), args );
  }
  cf_assert(m_faceDiffFlux.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::configureVolTermComputer( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFDMethodData> thisPtr(this);

  CFLog(INFO,"SpectralFD: Using volume term computer: " << m_volTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider(m_volTermComputerStr);
    cf_assert(prov.isNotNull());
    m_volTermComputer = prov->create(m_volTermComputerStr,thisPtr);
    configureNested ( m_volTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing BaseVolTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider("BaseVolTermComputer");
    cf_assert(prov.isNotNull());
    m_volTermComputer = prov->create("BaseVolTermComputer", thisPtr);
    configureNested ( m_volTermComputer.getPtr(), args );
  }
  cf_assert(m_volTermComputer.isNotNull());

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider(m_volTermComputerStr);
    cf_assert(prov.isNotNull());
    m_volTermComputer2nd = prov->create(m_volTermComputerStr,thisPtr);
    configureNested ( m_volTermComputer2nd.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing BaseVolTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider("BaseVolTermComputer");
    cf_assert(prov.isNotNull());
    m_volTermComputer2nd = prov->create("BaseVolTermComputer", thisPtr);
    configureNested ( m_volTermComputer2nd.getPtr(), args );
  }
  cf_assert(m_volTermComputer2nd.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::configureFaceTermComputer( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFDMethodData> thisPtr(this);

  CFLog(INFO,"SpectralFD: Using face term computer: " << m_faceTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseFaceTermComputer > > prov =
      Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider(m_faceTermComputerStr);
    cf_assert(prov.isNotNull());
    m_faceTermComputer = prov->create(m_faceTermComputerStr,thisPtr);
    configureNested ( m_faceTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing BaseFaceTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseFaceTermComputer > > prov =
      Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider("BaseFaceTermComputer");
    cf_assert(prov.isNotNull());
    m_faceTermComputer = prov->create("BaseFaceTermComputer", thisPtr);
    configureNested ( m_faceTermComputer.getPtr(), args );
  }
  cf_assert(m_faceTermComputer.isNotNull());

  // create additional face term computers
  /// @note for now, hard coded for the 3D case
  /// Works for 2D as well, but to many face term computers are created
  m_addFaceTermComputers.resize(10);
  m_addFaceTermComputersSP.resize(10);
  for (CFuint i = 0; i < 10; ++i)
  {
    try
    {
      Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseFaceTermComputer > > prov =
          Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider(m_faceTermComputerStr);
          cf_assert(prov.isNotNull());
      m_addFaceTermComputers[i] = prov->create(m_faceTermComputerStr,thisPtr);
      configureNested(m_addFaceTermComputers[i].getPtr(), args);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      CFLog(VERBOSE, "Choosing BaseFaceTermComputer instead ...\n");

      Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseFaceTermComputer > > prov =
          Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider("BaseFaceTermComputer");
          cf_assert(prov.isNotNull());
      m_addFaceTermComputers[i] = prov->create("BaseFaceTermComputer", thisPtr);
      configureNested(m_addFaceTermComputers[i].getPtr(), args);
    }
    cf_assert(m_addFaceTermComputers[i].isNotNull());
    m_addFaceTermComputersSP[i] = m_addFaceTermComputers[i].getPtr();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::configureBndFaceTermComputer( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFDMethodData> thisPtr(this);

  CFLog(INFO,"SpectralFD: Using boundary face term computer: " << m_bndFaceTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseBndFaceTermComputer > > prov =
      Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider(m_bndFaceTermComputerStr);
    cf_assert(prov.isNotNull());
    m_bndFaceTermComputer = prov->create(m_bndFaceTermComputerStr,thisPtr);
    configureNested ( m_bndFaceTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing BaseBndFaceTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseBndFaceTermComputer > > prov =
      Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider("BaseBndFaceTermComputer");
    cf_assert(prov.isNotNull());
    m_bndFaceTermComputer = prov->create("BaseBndFaceTermComputer", thisPtr);
    configureNested ( m_bndFaceTermComputer.getPtr(), args );
  }
  cf_assert(m_bndFaceTermComputer.isNotNull());

  // create additional face term computers
  /// @note for now, hard coded for the 3D case
  /// Works for 2D as well, but to much face term computers are created
  m_addBndFaceTermComputers.resize(10);
  m_addBndFaceTermComputersSP.resize(10);
  for (CFuint i = 0; i < 10; ++i)
  {
    try
    {
      Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseBndFaceTermComputer > > prov =
          Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider(m_bndFaceTermComputerStr);
          cf_assert(prov.isNotNull());
          m_addBndFaceTermComputers[i] = prov->create(m_bndFaceTermComputerStr,thisPtr);
          configureNested(m_addBndFaceTermComputers[i].getPtr(), args);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      CFLog(VERBOSE, "Choosing BaseBndFaceTermComputer instead ...\n");

      Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BaseBndFaceTermComputer > > prov =
          Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider("BaseBndFaceTermComputer");
          cf_assert(prov.isNotNull());
          m_addBndFaceTermComputers[i] = prov->create("BaseBndFaceTermComputer", thisPtr);
          configureNested(m_addBndFaceTermComputers[i].getPtr(), args);
    }
    cf_assert(m_addBndFaceTermComputers[i].isNotNull());
    m_addBndFaceTermComputersSP[i] = m_addBndFaceTermComputers[i].getPtr();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::configureBCStateComputers( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  cf_assert(m_bcTypeStr.size() == m_bcNameStr.size());

  SharedPtr<SpectralFDMethodData> thisPtr(this);

  // number of boundary conditions
  const CFuint nbrBcs = m_bcTypeStr.size();

  // resize m_bcs and m_bcsSP
  m_bcs.resize(nbrBcs);
  m_bcsSP.resize(nbrBcs);

  for (CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    CFLog(INFO, "BC type = " << m_bcTypeStr[iBc] << "\n");
    CFLog(INFO, "BC name = " << m_bcNameStr[iBc] << "\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFDMethodData , BCStateComputer > > prov =
      Environment::Factory<BCStateComputer>::getInstance().getProvider(m_bcTypeStr[iBc]);
    cf_assert(prov.isNotNull());
    m_bcs[iBc] = prov->create(m_bcNameStr[iBc],thisPtr);
    configureNested(m_bcs[iBc].getPtr(), args);
    cf_assert(m_bcs[iBc].isNotNull());

    // set SafePtr
    m_bcsSP[iBc] = m_bcs[iBc].getPtr();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::unsetup()
{
  SpaceMethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::setup()
{
  CFAUTOTRACE;

  SpaceMethodData::setup();

  // setup TRS Geo builder
  m_stdTrsGeoBuilder.setup();

  // setup face builder
  m_faceBuilder.setup();

  // setup cell builders
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

  // create local SV data
  createSDLocalData();

  // setup StatesReconstructor
  m_statesReconstructor->setup();

  // set the hasDiffTerm boolean
  /// @note it would be better to check a name related to the DiffusiveVarSet here
  m_hasDiffTerm = (_diffusiveVarStr != "Null");
  
  // create the transformer from update to solution variables ----------------
  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp =
  NamespaceSwitcher::getInstance().getNamespace(name);
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
  CFLog(INFO, "SpectralFDMethodData::setup() => updateToSolutionVarName = " << updateToSolutionVecTransStr << "\n");

  cf_assert(vecTransProv.isNotNull());
  m_updateToSolutionVecTrans.reset(vecTransProv->create(physModel->getImplementor()));
  cf_assert(m_updateToSolutionVecTrans.isNotNull());
  m_updateToSolutionVecTrans->setup(2);
  
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDMethodData::createSDLocalData()
{
  CFAUTOTRACE;

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // resize sdLocalData
  m_sdLocalData.resize(nbrElemTypes);

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
        throw Common::NotImplementedException (FromHere(),"Spectral difference has not been implemented for 1D");
      } break;
      case CFGeoShape::TRIAG:
      {
        throw Common::NotImplementedException (FromHere(),"Spectral difference has not been implemented for triangular cells");
      } break;
      case CFGeoShape::QUAD:
      {
        m_sdLocalData[iElemType] = new QuadSpectralFDElementData(polyOrder);
      } break;
      case CFGeoShape::TETRA:
      {
        throw Common::NotImplementedException (FromHere(),"Spectral difference has not been implemented for tetrahedral cells");
      } break;
      case CFGeoShape::HEXA:
      {
        m_sdLocalData[iElemType] = new HexaSpectralFDElementData(polyOrder);
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Spectral finite difference method not implemented for elements of type "
                                      + StringOps::to_str(elemShape) + ".");
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< ReconstructStatesSpectralFD > SpectralFDMethodData::getStatesReconstructor()
// this function has to be put in the implementation file because of the forward declaration of ReconstructStatesSpectralFD
{
  return m_statesReconstructor.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseBndFaceTermComputer > SpectralFDMethodData::getBndFaceTermComputer()
{
  return m_bndFaceTermComputer.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseBndFaceTermComputer >
    SpectralFDMethodData::getAdditionalBndFaceTermComputer(const CFuint idx)
{
  cf_assert(idx < m_addBndFaceTermComputers.size());
  return m_addBndFaceTermComputers[idx].getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< vector< SafePtr< BaseBndFaceTermComputer > > >
    SpectralFDMethodData::getAdditionalBndFaceTermComputers()
{
  return &m_addBndFaceTermComputersSP;
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseFaceTermComputer > SpectralFDMethodData::getFaceTermComputer()
{
  return m_faceTermComputer.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseFaceTermComputer >
    SpectralFDMethodData::getAdditionalFaceTermComputer(const CFuint idx)
{
  cf_assert(idx < m_addFaceTermComputers.size());
  return m_addFaceTermComputers[idx].getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< vector< SafePtr< BaseFaceTermComputer > > >
    SpectralFDMethodData::getAdditionalFaceTermComputers()
{
  return &m_addFaceTermComputersSP;
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseVolTermComputer > SpectralFDMethodData::getVolTermComputer()
{
  return m_volTermComputer.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseVolTermComputer > SpectralFDMethodData::getSecondVolTermComputer()
{
  return m_volTermComputer2nd.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< vector< SafePtr< BCStateComputer > > > SpectralFDMethodData::getBCStateComputers()
{
  return &m_bcsSP;
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< RiemannFlux > SpectralFDMethodData::getRiemannFlux()
{
  return m_riemannFlux.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

SafePtr< FaceDiffusiveFlux > SpectralFDMethodData::getFaceDiffusiveFlux()
{
  return m_faceDiffFlux.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SpectralFDElementData* >& SpectralFDMethodData::getSDLocalData()
{
  return m_sdLocalData;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD

