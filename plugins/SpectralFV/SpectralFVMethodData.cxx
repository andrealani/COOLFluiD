#include "Common/FilesystemException.hh"

#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SpectralFV/BaseBndFaceTermComputer.hh"
#include "SpectralFV/BaseFaceTermComputer.hh"
#include "SpectralFV/BaseVolTermComputer.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/FaceDiffusiveFlux.hh"
#include "SpectralFV/LineSpectralFVElementData.hh"
#include "SpectralFV/ReconstructStatesSpectralFV.hh"
#include "SpectralFV/RiemannFlux.hh"
#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SpectralFVElementData.hh"
#include "SpectralFV/SpectralFVMethodData.hh"
#include "SpectralFV/TetraSpectralFVElementData.hh"
#include "SpectralFV/TriagSpectralFVElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<
  NullMethodCommand< SpectralFVMethodData >,SpectralFVMethodData,SpectralFVModule >
  nullSpectralFVMethodComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::defineConfigOptions(Config::OptionList& options)
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
}

//////////////////////////////////////////////////////////////////////////////

SpectralFVMethodData::SpectralFVMethodData(Common::SafePtr<Framework::Method> owner) :
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
  m_svLocalData(),
  m_maxNbrStatesData(),
  m_maxNbrRFluxPnts(),
  m_hasDiffTerm(),
  m_resFactor(),
  m_createVolumesSocketBool()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_linearVarStr = "Roe";
  setParameter( "LinearVar", &m_linearVarStr );

  m_riemannFluxStr = "RoeFlux";
  setParameter("RiemannFlux", &m_riemannFluxStr);

  m_faceDiffFluxStr = "LocalApproach";
  setParameter("FaceDiffFlux", &m_faceDiffFluxStr);

  m_bndFaceTermComputerStr = "StdBndFaceTermComputer";
  setParameter("BndFaceTermComputer", &m_bndFaceTermComputerStr);

  m_faceTermComputerStr = "StdFaceTermComputer";
  setParameter("FaceTermComputer", &m_faceTermComputerStr);

  m_volTermComputerStr = "StdVolTermComputer";
  setParameter("VolTermComputer", &m_volTermComputerStr);

  m_createVolumesSocketBool = false;
  setParameter("ComputeVolumeForEachState", &m_createVolumesSocketBool);

  // options for bc commands
  m_bcTypeStr = vector<std::string>();
  setParameter("BcTypes",&m_bcTypeStr);

  m_bcNameStr = vector<std::string>();
  setParameter("BcNames",&m_bcNameStr);
}

//////////////////////////////////////////////////////////////////////////////

SpectralFVMethodData::~SpectralFVMethodData()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::configure ( Config::ConfigArgs& args )
{
  SpaceMethodData::configure(args);
  SharedPtr< SpectralFVMethodData > thisPtr(this);

  /* add here different strategies configuration */
  CFLog(INFO,"SpectralFVMethod: configureBndFaceTermComputer()\n");
  configureBndFaceTermComputer(args);

  CFLog(INFO,"SpectralFVMethod: configureFaceTermComputer()\n");
  configureFaceTermComputer(args);

  CFLog(INFO,"SpectralFVMethod: configureVolTermComputer()\n");
  configureVolTermComputer(args);

  CFLog(INFO,"SpectralFVMethod: configureBCStateComputers()\n");
  configureBCStateComputers(args);

  CFLog(INFO,"SpectralFVMethod: configureRiemannFlux()\n");
  configureRiemannFlux(args);

  CFLog(INFO,"SpectralFVMethod: configureFaceDiffusiveFlux()\n");
  configureFaceDiffusiveFlux(args);

  // states reconstructor
  Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , ReconstructStatesSpectralFV > > recProv =
    Environment::Factory< ReconstructStatesSpectralFV >::getInstance().getProvider("ReconstructStatesSpectralFV");
  cf_assert(recProv.isNotNull());
  m_statesReconstructor = recProv->create("ReconstructStatesSpectralFV",thisPtr);
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::configureRiemannFlux( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFVMethodData> thisPtr(this);

  CFLogDebugMin("SpectralFV: Using Riemann flux: " << m_riemannFluxStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , RiemannFlux > > prov =
      Environment::Factory< RiemannFlux >::getInstance().getProvider(m_riemannFluxStr);
    cf_assert(prov.isNotNull());
    m_riemannFlux = prov->create(m_riemannFluxStr,thisPtr);
    configureNested ( m_riemannFlux.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing RoeFlux instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , RiemannFlux > > prov =
      Environment::Factory< RiemannFlux >::getInstance().getProvider("RoeFlux");
    cf_assert(prov.isNotNull());
    m_riemannFlux = prov->create("RoeFlux", thisPtr);
    configureNested ( m_riemannFlux.getPtr(), args );
  }
  cf_assert(m_riemannFlux.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::configureFaceDiffusiveFlux( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFVMethodData> thisPtr(this);

  CFLogDebugMin("SpectralFV: Using face diffusive flux: " << m_faceDiffFluxStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , FaceDiffusiveFlux > > prov =
      Environment::Factory< FaceDiffusiveFlux >::getInstance().getProvider(m_faceDiffFluxStr);
    cf_assert(prov.isNotNull());
    m_faceDiffFlux = prov->create(m_faceDiffFluxStr,thisPtr);
    configureNested ( m_faceDiffFlux.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing RoeFlux instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , FaceDiffusiveFlux > > prov =
      Environment::Factory< FaceDiffusiveFlux >::getInstance().getProvider("LocalApproach");
    cf_assert(prov.isNotNull());
    m_faceDiffFlux = prov->create("LocalApproach", thisPtr);
    configureNested ( m_faceDiffFlux.getPtr(), args );
  }
  cf_assert(m_faceDiffFlux.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::configureVolTermComputer( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFVMethodData> thisPtr(this);

  CFLog(INFO,"SpectralFV: Using volume term computer: " << m_volTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider(m_volTermComputerStr);
    cf_assert(prov.isNotNull());
    m_volTermComputer = prov->create(m_volTermComputerStr,thisPtr);
    configureNested ( m_volTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing StdVolTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider("StdTermComputer");
    cf_assert(prov.isNotNull());
    m_volTermComputer = prov->create("StdVolTermComputer", thisPtr);
    configureNested ( m_volTermComputer.getPtr(), args );
  }
  cf_assert(m_volTermComputer.isNotNull());

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider(m_volTermComputerStr);
    cf_assert(prov.isNotNull());
    m_volTermComputer2nd = prov->create(m_volTermComputerStr,thisPtr);
    configureNested ( m_volTermComputer2nd.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing StdVolTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseVolTermComputer > > prov =
      Environment::Factory< BaseVolTermComputer >::getInstance().getProvider("StdTermComputer");
    cf_assert(prov.isNotNull());
    m_volTermComputer = prov->create("StdVolTermComputer", thisPtr);
    configureNested ( m_volTermComputer.getPtr(), args );
  }
  cf_assert(m_volTermComputer2nd.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::configureFaceTermComputer( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFVMethodData> thisPtr(this);

  CFLog(INFO,"SpectralFV: Using face term computer: " << m_faceTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseFaceTermComputer > > prov =
      Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider(m_faceTermComputerStr);
    cf_assert(prov.isNotNull());
    m_faceTermComputer = prov->create(m_faceTermComputerStr,thisPtr);
    configureNested ( m_faceTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing StdFaceTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseFaceTermComputer > > prov =
      Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider("StdFaceTermComputer");
    cf_assert(prov.isNotNull());
    m_faceTermComputer = prov->create("StdFaceTermComputer", thisPtr);
    configureNested ( m_faceTermComputer.getPtr(), args );
  }
  cf_assert(m_faceTermComputer.isNotNull());

  // create additional face term computers
  /// @note for now, hard coded for the 3D case
  /// Works for 2D as well, but to much face term computers are created
  m_addFaceTermComputers.resize(6);
  m_addFaceTermComputersSP.resize(6);
  for (CFuint i = 0; i < 6; ++i)
  {
    try
    {
      Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseFaceTermComputer > > prov =
          Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider(m_faceTermComputerStr);
          cf_assert(prov.isNotNull());
      m_addFaceTermComputers[i] = prov->create(m_faceTermComputerStr,thisPtr);
      configureNested(m_addFaceTermComputers[i].getPtr(), args);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      CFLog(VERBOSE, "Choosing StdFaceTermComputer instead ...\n");

      Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseFaceTermComputer > > prov =
          Environment::Factory<BaseFaceTermComputer>::getInstance().getProvider("StdFaceTermComputer");
          cf_assert(prov.isNotNull());
      m_addFaceTermComputers[i] = prov->create("StdFaceTermComputer", thisPtr);
      configureNested(m_addFaceTermComputers[i].getPtr(), args);
    }
    cf_assert(m_addFaceTermComputers[i].isNotNull());
    m_addFaceTermComputersSP[i] = m_addFaceTermComputers[i].getPtr();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::configureBndFaceTermComputer ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SpectralFVMethodData> thisPtr(this);

  CFLog(INFO,"SpectralFV: Using boundary face term computer: " << m_bndFaceTermComputerStr << "\n");

  try
  {
    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseBndFaceTermComputer > > prov =
      Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider(m_bndFaceTermComputerStr);
    cf_assert(prov.isNotNull());
    m_bndFaceTermComputer = prov->create(m_bndFaceTermComputerStr,thisPtr);
    configureNested ( m_bndFaceTermComputer.getPtr(), args );
  }
  catch (Common::NoSuchValueException& e)
  {
    CFLog(VERBOSE, e.what() << "\n");
    CFLog(VERBOSE, "Choosing StdBndFaceTermComputer instead ...\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseBndFaceTermComputer > > prov =
      Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider("StdBndFaceTermComputer");
    cf_assert(prov.isNotNull());
    m_bndFaceTermComputer = prov->create("StdBndFaceTermComputer", thisPtr);
    configureNested ( m_bndFaceTermComputer.getPtr(), args );
  }
  cf_assert(m_bndFaceTermComputer.isNotNull());

  // create additional face term computers
  /// @note for now, hard coded for the 3D case
  /// Works for 2D as well, but to much face term computers are created
  m_addBndFaceTermComputers.resize(6);
  m_addBndFaceTermComputersSP.resize(6);
  for (CFuint i = 0; i < 6; ++i)
  {
    try
    {
      Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseBndFaceTermComputer > > prov =
          Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider(m_bndFaceTermComputerStr);
          cf_assert(prov.isNotNull());
          m_addBndFaceTermComputers[i] = prov->create(m_bndFaceTermComputerStr,thisPtr);
          configureNested(m_addBndFaceTermComputers[i].getPtr(), args);
    }
    catch (Common::NoSuchValueException& e)
    {
      CFLog(VERBOSE, e.what() << "\n");
      CFLog(VERBOSE, "Choosing StdBndFaceTermComputer instead ...\n");

      Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BaseBndFaceTermComputer > > prov =
          Environment::Factory<BaseBndFaceTermComputer>::getInstance().getProvider("StdBndFaceTermComputer");
          cf_assert(prov.isNotNull());
          m_addBndFaceTermComputers[i] = prov->create("StdBndFaceTermComputer", thisPtr);
          configureNested(m_addBndFaceTermComputers[i].getPtr(), args);
    }
    cf_assert(m_addBndFaceTermComputers[i].isNotNull());
    m_addBndFaceTermComputersSP[i] = m_addBndFaceTermComputers[i].getPtr();
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::configureBCStateComputers( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  cf_assert(m_bcTypeStr.size() == m_bcNameStr.size());

  SharedPtr<SpectralFVMethodData> thisPtr(this);

  // number of boundary conditions
  const CFuint nbrBcs = m_bcTypeStr.size();

  // resize m_bcs and m_bcsSP
  m_bcs.resize(nbrBcs);
  m_bcsSP.resize(nbrBcs);

  for (CFuint iBc = 0; iBc < nbrBcs; ++iBc)
  {
    CFLog(INFO, "BC type = " << m_bcTypeStr[iBc] << "\n");
    CFLog(INFO, "BC name = " << m_bcNameStr[iBc] << "\n");

    Common::SafePtr<BaseMethodStrategyProvider< SpectralFVMethodData , BCStateComputer > > prov =
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

void SpectralFVMethodData::unsetup()
{
  SpaceMethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::setup()
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
  /// @note think about making the creation of the NumericalJacobian object optional
  /// (only if the Jacobian has to be computed)
  m_numJacob.reset(new NumericalJacobian("NumericalJacobian"));

  // set reference values in numerical Jacobian computer
  RealVector refValues = PhysicalModelStack::getActive()->getImplementor()->getRefStateValues();
  m_numJacob->setRefValues(refValues);

  // setup the variable sets
  _updateVar   ->setup();
  _solutionVar ->setup();
  _diffusiveVar->setup();

  // create local SV data
  createSVLocalData();

  // setup StatesReconstructor
  m_statesReconstructor->setup();

  // set the hasDiffTerm boolean
  /// @note it would be better to check a name related to the DiffusiveVarSet here
  m_hasDiffTerm = (_diffusiveVarStr != "Null");
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVMethodData::createSVLocalData()
{
  CFAUTOTRACE;

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // resize svLocalData
  m_svLocalData.resize(nbrElemTypes);

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
        m_svLocalData[iElemType] = new LineSpectralFVElementData(polyOrder);
      } break;
      case CFGeoShape::TRIAG:
      {
        m_svLocalData[iElemType] = new TriagSpectralFVElementData(polyOrder);
      } break;
      case CFGeoShape::TETRA:
      {
        m_svLocalData[iElemType] = new TetraSpectralFVElementData(polyOrder);
      } break;
      default:
      {
        throw Common::NotImplementedException (FromHere(),"Spectral finite volume method only implemented for elements of type LINE (1D), TRIAG (2D) and TETRA (3D).");
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< ReconstructStatesSpectralFV > SpectralFVMethodData::getStatesReconstructor()
// this function has to be put in the implementation file because of the forward declaration of ReconstructStatesSpectralFV
{
  return m_statesReconstructor.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseBndFaceTermComputer > SpectralFVMethodData::getBndFaceTermComputer()
{
  return m_bndFaceTermComputer.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseBndFaceTermComputer >
    SpectralFVMethodData::getAdditionalBndFaceTermComputer(const CFuint idx)
{
  cf_assert(idx < m_addBndFaceTermComputers.size());
  return m_addBndFaceTermComputers[idx].getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< vector< SafePtr< BaseBndFaceTermComputer > > >
    SpectralFVMethodData::getAdditionalBndFaceTermComputers()
{
  return &m_addBndFaceTermComputersSP;
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseFaceTermComputer > SpectralFVMethodData::getFaceTermComputer()
{
  return m_faceTermComputer.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseFaceTermComputer >
    SpectralFVMethodData::getAdditionalFaceTermComputer(const CFuint idx)
{
  cf_assert(idx < m_addFaceTermComputers.size());
  return m_addFaceTermComputers[idx].getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< vector< SafePtr< BaseFaceTermComputer > > >
    SpectralFVMethodData::getAdditionalFaceTermComputers()
{
  return &m_addFaceTermComputersSP;
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseVolTermComputer > SpectralFVMethodData::getVolTermComputer()
{
  return m_volTermComputer.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< BaseVolTermComputer > SpectralFVMethodData::getSecondVolTermComputer()
{
  return m_volTermComputer2nd.getPtr();
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< vector< SafePtr< BCStateComputer > > > SpectralFVMethodData::getBCStateComputers()
{
  return &m_bcsSP;
}

/////////////////////////////////////////////////////////////////////////////

SafePtr< RiemannFlux > SpectralFVMethodData::getRiemannFlux()
{
  return m_riemannFlux.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

SafePtr< FaceDiffusiveFlux > SpectralFVMethodData::getFaceDiffusiveFlux()
{
  return m_faceDiffFlux.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD

