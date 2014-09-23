
#include "SGSViscosity.hh"
#include "Framework/DataHandle.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "LESDataProcessing.hh"
#include "Framework/ConvectiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<SGSViscosity, 
                       LESProcessingData, 
                       TurbulenceFunction, 
                       LESDataProcessingModule> 
sgsViscosityProvider("SGSViscosity");

//////////////////////////////////////////////////////////////////////////////

void SGSViscosity::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

SGSViscosity::SGSViscosity(const std::string& name) :
  TurbulenceFunction(name)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.  
}

//////////////////////////////////////////////////////////////////////////////


SGSViscosity::~SGSViscosity()
{
}

//////////////////////////////////////////////////////////////////////////////

void SGSViscosity::configure ( Config::ConfigArgs& args )
{
  TurbulenceFunction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SGSViscosity::setup()
{
  CFAUTOTRACE; 
  TurbulenceFunction::setup();

  // Set the size of the socket
  m_socket.resize(getNbStates());
  
  m_nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  m_gradientPtrs.resize(m_nbEqs);
  m_dimState.resize(m_nbEqs);  
  m_lesVar = getMethodData().getLESVar().d_castTo<LESVAR>();
  
}

//////////////////////////////////////////////////////////////////////////////

void SGSViscosity::compute(const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell)
{
  for(CFuint iEq=0; iEq<m_nbEqs; iEq++)   m_gradientPtrs[iEq] = &gradients[iEq];
  getMethodData().getUpdateVarSet()->setDimensionalValues(state, m_dimState);
  m_socket[iCell] = m_lesVar->getDynSGSViscosity(m_dimState,m_gradientPtrs,getMethodData().getVolume(iCell));  
}
      
//////////////////////////////////////////////////////////////////////////////

void SGSViscosity::computeAverage(const CFreal& oldWeight, const CFreal& newWeight, const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell)
{
  for(CFuint iEq=0; iEq<m_nbEqs; iEq++)   m_gradientPtrs[iEq] = &gradients[iEq];
  getMethodData().getUpdateVarSet()->setDimensionalValues(state, m_dimState);
  m_socket[iCell] = oldWeight*m_socket[iCell] + newWeight*m_lesVar->getDynSGSViscosity(m_dimState,m_gradientPtrs,getMethodData().getVolume(iCell));  
}

//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
