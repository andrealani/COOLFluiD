#include "LESDataProcessing/LESDataProcessing.hh"
#include "ComputeTurbulenceFunctions.hh"
#include "Framework/MethodCommandProvider.hh"
#include "TurbulenceFunction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeTurbulenceFunctions, 
                      LESProcessingData, 
                      LESDataProcessingModule> 
LESDataProcessingComputeTurbulenceFunctionsProvider("ComputeTurbulenceFunctions");

//////////////////////////////////////////////////////////////////////////////

void ComputeTurbulenceFunctions::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("ProcessRate","Rate to process the data.");
}

//////////////////////////////////////////////////////////////////////////////

ComputeTurbulenceFunctions::ComputeTurbulenceFunctions(const std::string& name) 
: LESProcessingCom(name)
{
  addConfigOptionsTo(this);
  // by default the data processing is run once -> processRate=infinity
  m_processRate = MathTools::MathConsts::CFuintMax();
  setParameter("ProcessRate",&m_processRate);
}

////////////////////////////////////////////////////////////////////////////////

void ComputeTurbulenceFunctions::configure ( Config::ConfigArgs& args )
{
  LESProcessingCom::configure(args);
  m_globalSockets.createSocketSink<Framework::State*>("states");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeTurbulenceFunctions::setup()
{
  CFLog(INFO, " +++ ComputeTurbulenceFunctions::setup() \n");
  DataHandle<Framework::State*,Framework::GLOBAL> states = m_globalSockets.getSocketSink<Framework::State*>("states")->getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
    
  m_gradients.resize(nbEqs);
  for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
    m_gradients[iEq].resize(dim);
  }
  
  m_turbulenceFunctions = getMethodData().getTurbulenceFunctions();
  
}

//////////////////////////////////////////////////////////////////////////////

void ComputeTurbulenceFunctions::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeTurbulenceFunctions::execute()
{
  
  if((!(Framework::SubSystemStatusStack::getActive()->getNbIter() % m_processRate)) ||
    (Framework::SubSystemStatusStack::getActive()->getNbIter() == 1)) {
      
    CFLog(INFO, "Calculating Turbulence Structures \n");
    
    DataHandle<Framework::State*,Framework::GLOBAL> states = m_globalSockets.getSocketSink<Framework::State*>("states")->getDataHandle();
    const CFuint nbCells = states.size();       
    const CFuint nbTurbulenceFunctions = m_turbulenceFunctions.size();
    
    // Calculate gradients for every cell and fill the source sockets
    for(CFuint iCell=0; iCell<nbCells; ++iCell) {
      
      getMethodData().getGradientComputer()->compute(m_gradients,iCell);
      
      for(CFuint iFunc=0; iFunc<nbTurbulenceFunctions; ++iFunc) {
        m_turbulenceFunctions[iFunc]->compute((*states[iCell]),m_gradients,iCell);
      }

    }
  }
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESDataProcessing

  } // namespace Numerics

} // namespace COOLFluiD
