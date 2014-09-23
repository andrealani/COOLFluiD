
#include "CFL.hh"
#include "Framework/DataHandle.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "LESDataProcessing.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MeshData.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<CFL, 
                       LESProcessingData, 
                       TurbulenceFunction, 
                       LESDataProcessingModule> 
CFLsocketProvider("CFL");

//////////////////////////////////////////////////////////////////////////////

void CFL::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("DT","Time step to calculate CFL with. (default = 1.0)");
}

//////////////////////////////////////////////////////////////////////////////

CFL::CFL(const std::string& name) :
  TurbulenceFunction(name),
  m_updateCoeff(CFNULL)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.  
  m_dt = 1.0;
  setParameter("DT",&m_dt);
}

//////////////////////////////////////////////////////////////////////////////


CFL::~CFL()
{
}

//////////////////////////////////////////////////////////////////////////////

void CFL::configure ( Config::ConfigArgs& args )
{
  TurbulenceFunction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CFL::setup()
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
  
  m_dt /= (PhysicalModelStack::getActive()->getImplementor()->getRefTime());
  
  
}

//////////////////////////////////////////////////////////////////////////////

void CFL::compute(const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell)
{
  
  // Calculate DIMENSIONAL primitive state
  // m_primState = getMethodData().transformToPrimDim(state);
  
  CFreal volume = getMethodData().getVolumeAdim(iCell);
  // m_socket[iCell] = m_updateCoeff[iCell] ;
  
  m_socket[iCell] = m_dt * m_updateCoeff[iCell] / volume;
}

//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
