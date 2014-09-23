#include "STM_LDASchemeSysHO.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "MathTools/MatrixInverter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STM_LDASchemeSysHO,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeLDASchemeSyshoProvider("STM_SysLDAC_HO");

//////////////////////////////////////////////////////////////////////////////

STM_LDASchemeSysHO::STM_LDASchemeSysHO(const std::string& name) :
  STM_SplitterSys(name),
  _sumKplus(),
  _invK(),
  m_temp_mat(),
m_utemp()
{
}

//////////////////////////////////////////////////////////////////////////////

STM_LDASchemeSysHO::~STM_LDASchemeSysHO()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeSysHO::setup()
{
  STM_SplitterSys::setup();

  _sumKplus.resize(_nbEquations,_nbEquations);
  _invK.resize(_nbEquations,_nbEquations);
  m_temp_mat.resize(_nbEquations,_nbEquations);
m_utemp.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeSysHO::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"STM_LDASchemeScalarHO::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////
void STM_LDASchemeSysHO::distributePast(vector<RealVector>& residual )
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
 
  DistributionData& ddata = getMethodData().getDistributionData();

  RealVector m_phi = ddata.phi;

  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKplus  += *_kPlus[iState];
  }
 
  _inverter->invert(_sumKplus, _invK);

  vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
  m_temp_mat = *_kPlus[iState]*_invK;

m_temp_mat *= timeStep;

m_utemp = m_temp_mat*m_phi;

      residual[iState] = m_utemp;


      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = *_kPlus[iState]*_invK;

    }
}


//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeSysHO::distribute(vector<RealVector>& residual)
{
  // Compute sumKmin
  
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();

 DistributionData& ddata = getMethodData().getDistributionData();

  RealVector m_phi = ddata.phi;

  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKplus  += *_kPlus[iState];
  }
 

_inverter->invert(_sumKplus, _invK);

  vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
  m_temp_mat = *_kPlus[iState]*_invK;

m_temp_mat *= timeStep;

m_utemp = m_temp_mat*m_phi;

      residual[iState] = m_utemp;

      betas[iState] = *_kPlus[iState]*_invK;

     

    }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
