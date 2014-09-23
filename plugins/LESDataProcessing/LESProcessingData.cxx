#include "LESDataProcessing.hh"
#include "LESProcessingData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
// #include "Framework/PhysicalModelImpl.hh"
// #include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "LES/LESVarSet.hh"
#include "GradientComputer.hh"
#include "TurbulenceFunction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<LESProcessingData>,
											LESProcessingData,
											LESDataProcessingModule>
nullLESProcessingComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void LESProcessingData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("GradientComputer","Which gradient computer should be used (default=GradientComputerFVMCC)");
  options.addConfigOption< std::vector<std::string> >("TurbulenceFunctions","Which Turbulence Functions are to be used.");
  options.addConfigOption< std::vector<std::string> >("AverageTurbulenceFunctions","Which Turbulence Functions are to be used to average.");
}

//////////////////////////////////////////////////////////////////////////////

LESProcessingData::LESProcessingData(Common::SafePtr<Framework::Method> owner)
  : DataProcessingData(owner),
    m_gradientComputer(),
    m_turbulenceFunctions(),
    m_averageTurbulenceFunctions()
{
   addConfigOptionsTo(this);
   m_gradientComputerStr = "GradientComputerFVMCC";
   setParameter("GradientComputer",&m_gradientComputerStr);
   
   m_turbulenceFunctionNames = std::vector<std::string>();
   setParameter("TurbulenceFunctions",&m_turbulenceFunctionNames);
   
   m_averageTurbulenceFunctionNames = std::vector<std::string>();
   setParameter("AverageTurbulenceFunctions",&m_averageTurbulenceFunctionNames);
}

//////////////////////////////////////////////////////////////////////////////

LESProcessingData::~LESProcessingData()
{
}

//////////////////////////////////////////////////////////////////////////////

void LESProcessingData::configure ( Config::ConfigArgs& args )
{
  DataProcessingData::configure(args);
  
  
  m_updateVar = getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVar();
  
  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getNamespace();
  Common::SafePtr<Namespace> nsp = Framework::NamespaceSwitcher::getInstance().getNamespace(namespc);
  Common::SafePtr<PhysicalModel> physModel = Framework::PhysicalModelStack::getInstance().getEntryByNamespace(nsp);


  std::string updateVarStr = getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVarStr();

  // create the transformer from update to primitive variables
  std::string provider =
      Framework::VarSetTransformer::getProviderName(physModel->getConvectiveName(), updateVarStr, "Prim");
      
  CFLog(VERBOSE, "   +++ Creating VarSetTransformer: " << provider << " \n");
  m_updateToPrimVar =
      Environment::Factory<Framework::VarSetTransformer>::getInstance()
      .getProvider(provider)->create(physModel->getImplementor());
  cf_assert(m_updateToPrimVar.isNotNull());
  
  
  Common::SharedPtr<LESProcessingData> thisPtr(this);
    
  CFLog(INFO, " +++ Strategy " << m_gradientComputerStr << "::configure()\n");
  configureStrategy(args,m_gradientComputer,m_gradientComputerStr,m_gradientComputerStr,thisPtr);
  
  const CFuint nbTurbulenceFunctions = m_turbulenceFunctionNames.size();
  m_turbulenceFunctions.resize(nbTurbulenceFunctions);
  for(CFuint i=0; i<nbTurbulenceFunctions; ++i) {
    CFLog(INFO, " +++ TurbulenceFunction Strategy " << m_turbulenceFunctionNames[i] << "::configure()\n");
    configureStrategy(args,
                      m_turbulenceFunctions[i],
                      m_turbulenceFunctionNames[i],
                      m_turbulenceFunctionNames[i],
                      thisPtr);
  }
  
  const CFuint nbAverageTurbulenceFunctions = m_averageTurbulenceFunctionNames.size();
  m_averageTurbulenceFunctions.resize(nbAverageTurbulenceFunctions);
  for(CFuint i=0; i<nbAverageTurbulenceFunctions; ++i) {
    CFLog(INFO, " +++ Average TurbulenceFunction Strategy " << m_averageTurbulenceFunctionNames[i] << "::configure()\n");
    configureStrategy(args,
                      m_averageTurbulenceFunctions[i],
                      m_averageTurbulenceFunctionNames[i],
                      m_averageTurbulenceFunctionNames[i],
                      thisPtr);
  }

  
}

//////////////////////////////////////////////////////////////////////////////

void LESProcessingData::setup()
{
  DataProcessingData::setup();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // m_varSet = getUpdateVarSet().d_castTo<EULERVAR>();
  // setup the update to primitive variables transformer
  m_dimState.resize(nbEqs);  
  m_primState = new Framework::State(RealVector(nbEqs));
  m_updateToPrimVar->setup(1);
  // m_lesVarSet->setup();
    
}

//////////////////////////////////////////////////////////////////////////////

void LESProcessingData::unsetup()
{
  deletePtr(m_primState);
  DataProcessingData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector LESProcessingData::transformToPrim(RealVector& state)
{
  m_dimState = *getUpdateToPrimTransformer()->transform(static_cast<Framework::State*>(&state));
  return m_dimState;
}

//////////////////////////////////////////////////////////////////////////////

RealVector LESProcessingData::transformToPrimDim(const RealVector& state)
{
  m_updateVar->setDimensionalValues(state, *m_primState);
  m_dimState = *getUpdateToPrimTransformer()->transform(m_primState);
  return m_dimState;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::DiffusiveVarSet> LESProcessingData::getLESVar()
{
  return getCollaborator<SpaceMethod>()->getSpaceMethodData()->getDiffusiveVar();
}

//////////////////////////////////////////////////////////////////////////////

CFreal LESProcessingData::getVolume(const CFuint& cellID) {
  return m_gradientComputer->getVolume(cellID);
}

//////////////////////////////////////////////////////////////////////////////

CFreal LESProcessingData::getVolumeAdim(const CFuint& cellID) {
  return m_gradientComputer->getVolumeAdim(cellID);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESProcessing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

