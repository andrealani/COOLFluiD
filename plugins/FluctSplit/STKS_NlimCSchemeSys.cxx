#include "STKS_NlimCSchemeSys.hh"
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

MethodStrategyProvider<STKS_NlimCSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
stksNlimcSchemeSysProvider("STKS_SysNlimC");

////////////////////////////////////////////////////////////////////////////

void STKS_NlimCSchemeSys::defineConfigOptions(Config::OptionList& options)
{
  ///@todo NV: for 3D we should put as input two angles
   options.addConfigOption<CFreal>("ProjDir2D","Choosea direction to project the residual(an angle in degree)");
}
//////////////////////////////////////////////////////////////////////////////

STKS_NlimCSchemeSys::STKS_NlimCSchemeSys(const std::string& name) :
  STKS_SplitterSys(name),
  m_sumKplus(),
  _invK(),
  m_uTemp(),
  m_diss(),
  m_sumKplusU()
  { 
  addConfigOptionsTo(this);

  m_angle = 10.0;
  setParameter("ProjDir2D",&m_angle);
}

//////////////////////////////////////////////////////////////////////////////

STKS_NlimCSchemeSys::~STKS_NlimCSchemeSys()
{
}


//////////////////////////////////////////////////////////////////////////////

void STKS_NlimCSchemeSys::setup()
{
  STKS_SplitterSys::setup();

  getMethodData().getDistributionData().needDiss = true;
  
  m_sumKplus.resize(_nbEquations, _nbEquations);
  m_uTemp.resize(_nbEquations);
  m_sumKplusU.resize(_nbEquations);
  m_diss.resize(_nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _phy.resize(_nbEquations);
  _phyChar.resize(_nbEquations);
  _rightEigenVector.resize(_nbEquations, _nbEquations);
  _leftEigenVector.resize(_nbEquations, _nbEquations);
  _sumBeta.resize(_nbEquations);

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  past_diss.resize(maxNbStatesInCell);
  _residualChar.resize(maxNbStatesInCell);
  _beta.resize(maxNbStatesInCell);
  _betaLim.resize(maxNbStatesInCell);
  m_time_comp.resize(maxNbStatesInCell);
  
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){
   
    past_diss[iState].resize(_nbEquations);
    m_time_comp[iState].resize(_nbEquations);
    _residualChar[iState].resize(_nbEquations);
    _beta[iState].resize(_nbEquations);
    _betaLim[iState].resize(_nbEquations);
    m_time_comp[iState].resize(_nbEquations);
  }
  
  // we size the direction to 2 because 3D is not implemented for this
  m_normal.resize(2);

  m_normal[XX] = cos(m_angle*3.1415/180.);
  m_normal[YY] = sin(m_angle*3.1415/180.);

}

//////////////////////////////////////////////////////////////////////////////

void STKS_NlimCSchemeSys::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}


//////////////////////////////////////////////////////////////////////////////
  void STKS_NlimCSchemeSys::ComputePastDissipationAndTimeComp( const vector<State*>& tStates )
{
  const CFuint nbStatesInCell = _nbStatesInCell;

  m_sumKplus = *_kPlus[0];
  m_sumKplusU = *_kPlus[0]*(*tStates[0]);
  
  for (CFuint iState = 1; iState < nbStatesInCell; ++ iState){
    m_sumKplus += *_kPlus[iState];
    m_sumKplusU += *_kPlus[iState]*(*tStates[iState]);
  } 

  _inverter->invert(m_sumKplus, _invK);
  m_uTemp = _invK*m_sumKplusU;

  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    past_diss[iState]  = *_kPlus[iState]*(*tStates[iState])- (*_kPlus[iState])*m_uTemp;
    m_time_comp[iState] = (-1.)*(*tStates[iState]);
  }
}


//////////////////////////////////////////////////////////////////////////////

void STKS_NlimCSchemeSys::distribute(vector<RealVector>& residual)
{ 
  DistributionData& ddata = getMethodData().getDistributionData();
  const RealVector& m_phiT = ddata.phi;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const vector<State*>& tStates = *ddata.tStates;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT()/dim;
  const CFreal Area = _cellVolume/(dim+1.);

  m_sumKplus = *_kPlus[0]; 
  m_sumKplusU = *_kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    m_sumKplus  += *_kPlus[iState];
    m_sumKplusU += *_kPlus[iState]*(*tStates[iState]);
  }

  _inverter->invert(m_sumKplus, _invK);
  m_uTemp = _invK*m_phiT;
  
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    residual[iState] = *_kPlus[iState]*m_uTemp;
  }

  m_uTemp = _invK*m_sumKplusU;
 
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    // And then the dissipation
    m_diss = *_kPlus[iState]*(*tStates[iState]) - *_kPlus[iState]*m_uTemp;

   residual[iState] += timeStep*(m_diss + past_diss[iState]);
  
  residual[iState] += Area*((*tStates[iState]) +  m_time_comp[iState]);
  }

  m_normal[XX] = cos(m_angle*3.1415/180.);
  m_normal[YY] = sin(m_angle*3.1415/180.);

 // Compute the cell residual...
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _phy[iEq] = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      _phy[iEq] += residual[iState][iEq];
    }
  }

  // Transform into characteristic variables the residual, phy
  // Compute the EigenVectors
  getMethodData().getDistribVar()->computeEigenValuesVectors(_rightEigenVector,
							 _leftEigenVector,
							 _eValues,
							 m_normal);

  //Transform the residual into characteristic variables
  _phyChar = _leftEigenVector *_phy;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

      _residualChar[iState] = _leftEigenVector * residual[iState];

  }
  // Do the limiting just as in scalar case
  
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    if (std::abs(_phyChar[iEq])>MathTools::MathConsts::CFrealEps()){
      _sumBeta[iEq] = 0.;
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

        _beta[iState][iEq] = _residualChar[iState][iEq]/_phyChar[iEq];

        _sumBeta[iEq] += max(0.,_beta[iState][iEq]);

      }
    }
  }
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    if (std::abs(_phyChar[iEq])>MathTools::MathConsts::CFrealEps()){
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _betaLim[iState][iEq] = max(0.,_beta[iState][iEq])/_sumBeta[iEq];
        _residualChar[iState][iEq] = _betaLim[iState][iEq] * _phyChar[iEq];
        }
    }
    else{
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _residualChar[iState][iEq] = 0.;
      }
    }
  }

  //Transform the residual back into conservative variables
 _phy = _rightEigenVector * _phyChar;
 for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
     residual[iState] = _rightEigenVector * _residualChar[iState];
 }
  
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations

      RealMatrix& betaLDA = (*getMethodData().getDistributionData().currBetaMat)[iState];
      betaLDA = (*_kPlus[iState])*_invK;
    }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
