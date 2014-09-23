#include "STM_LDASchemeScalarHO.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STM_LDASchemeScalarHO,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeLDASchemeScalarhoProvider("STM_ScalarLDAC_HO");

//////////////////////////////////////////////////////////////////////////////

STM_LDASchemeScalarHO::STM_LDASchemeScalarHO(const std::string& name) :
  STM_SplitterScalar(name),
  _sumKplus()
{
}

//////////////////////////////////////////////////////////////////////////////

STM_LDASchemeScalarHO::~STM_LDASchemeScalarHO()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeScalarHO::setup()
{
  STM_SplitterScalar::setup();

  _sumKplus.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeScalarHO::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"STM_LDASchemeScalarHO::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////
void STM_LDASchemeScalarHO::distributePast(vector<RealVector>& residual )
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
 
  DistributionData& ddata = getMethodData().getDistributionData();

  RealVector m_phi = ddata.phi;

  _sumKplus = _kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKplus  += _kPlus[iState];
  }
 



  vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = _kPlus[iState]/_sumKplus;
      residual[iState] = (timeStep)*betas[iState]*m_phi;

    }
}


//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeScalarHO::distribute(vector<RealVector>& residual)
{
  // Compute sumKmin
  
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();

 DistributionData& ddata = getMethodData().getDistributionData();

  RealVector m_phi = ddata.phi;

  _sumKplus = _kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKplus  += _kPlus[iState];
  }
 



  vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = _kPlus[iState]/_sumKplus;

      residual[iState] = (timeStep)*betas[iState]*m_phi;
     

    }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
