
#include "LocalTimeStep.hh"
#include "Framework/DataHandle.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "LESDataProcessing.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MeshData.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/CFL.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<LocalTimeStep, 
                       LESProcessingData, 
                       TurbulenceFunction, 
                       LESDataProcessingModule> 
LocalTimeStepSocketProvider("LocalTimeStep");

//////////////////////////////////////////////////////////////////////////////

void LocalTimeStep::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

LocalTimeStep::LocalTimeStep(const std::string& name) :
  TurbulenceFunction(name),
  m_updateCoeff(CFNULL)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.  
}

//////////////////////////////////////////////////////////////////////////////


LocalTimeStep::~LocalTimeStep()
{
}

//////////////////////////////////////////////////////////////////////////////

void LocalTimeStep::configure ( Config::ConfigArgs& args )
{
  TurbulenceFunction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void LocalTimeStep::setup()
{
  CFAUTOTRACE; 
  TurbulenceFunction::setup();

  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  std::string datahandleName = nsp + "_updateCoeff";
  m_updateCoeff =
    MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

  // Set the size of the socket
  m_socket.resize(getNbStates());  
   
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  m_primState.resize(nbEqs);    
  
}

//////////////////////////////////////////////////////////////////////////////

void LocalTimeStep::compute(const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell)
{
  
  // Calculate DIMENSIONAL primitive state
  // m_primState = getMethodData().transformToPrimDim(state);
  
  const CFreal cfl = getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->getCFL()->getCFLValue();
  CFreal volume = getMethodData().getVolumeAdim(iCell);
  // m_socket[iCell] = m_updateCoeff[iCell] ;
  
  m_socket[iCell] = volume * cfl / m_updateCoeff[iCell] ;
}

//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
