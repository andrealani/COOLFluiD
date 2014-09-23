#include "STKT_SplitterScalarQuad.hh"
#include "Common/CFLog.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

STKT_SplitterScalarQuad::STKT_SplitterScalarQuad(const std::string& name) :
  SpaceTime_Splitter(name),
  _kPlus(0),
  _kMin(0),
  _k(0),
  _kPlus_space(0)
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_SplitterScalarQuad::~STKT_SplitterScalarQuad()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_SplitterScalarQuad::setBlockData()
{
  const CFuint totalNbEq = PhysicalModelStack::getActive()->getNbEq();

  cf_assert(_blockSeparator <= totalNbEq);

  if (_blockSeparator == totalNbEq) {
    _firstVarID = 0;
    _lastVarID  = totalNbEq;
    _nbEquations = totalNbEq;
  }
  else {
    _firstVarID = _blockSeparator;
    _lastVarID = totalNbEq;
    _nbEquations = _lastVarID - _firstVarID;
  }
  cf_assert(_nbEquations == _lastVarID - _firstVarID);
  cf_assert(_lastVarID <= totalNbEq);
}

//////////////////////////////////////////////////////////////////////////////

void STKT_SplitterScalarQuad::setup()
{
  SpaceTime_Splitter::setup();

  _isOnlySplitter = (_nbEquations == PhysicalModelStack::getActive()->getNbEq());

  CFLogDebugMax( "SpaceTimeRDS_SplitterScalar::_nbEquations: " <<
        _nbEquations << ", _firstVar = " << _firstVarID  <<
        ", _lastVar = " << _lastVarID << "\n");

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _kPlus.resize(maxNbStatesInCell);
  _kMin.resize(maxNbStatesInCell);
  _k.resize(maxNbStatesInCell);
  _kPlus_space.resize(maxNbStatesInCell);

  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _kPlus[i].resize(_nbEquations);
    _kMin[i].resize(_nbEquations);
    _k[i].resize(_nbEquations);
    _kPlus_space[i].resize(_nbEquations);
  }

  DIM = PhysicalModelStack::getActive()->getDim();
}

//////////////////////////////////////////////////////////////////////////////

void STKT_SplitterScalarQuad::computeK(const vector<State*>& states,
					   const InwardNormalsData* const normalsData)
{

  _nbStatesInCell = states.size();
  m_normals = normalsData;
  _nodeArea.resize(states.size());

   const CFreal kCoeff = 0.5;
  CFreal Area = _cellVolume;
  CFreal PastArea = _pastCellVolume;
 const CFreal tCoeff = _timeStep/2.;
  CFreal mini = 10.e+10;
  CFreal temp;
  CFreal oneforthArea = _cellVolume/4.0;



  //The transformation of the normal is needed if the coordinate system is rotated
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    SpaceTime_Splitter::setAdimensionalNormal(iState);

    getMethodData().getDistribVar()->computeScalarJacobian(_adimNormal, _k[iState]);

    _nodeArea[iState] = m_normals->getAreaNode(iState);
    _k[iState] *= _nodeArea[iState] * kCoeff;

    for (CFuint iEq = 0; iEq < _nbEquations; iEq++) {
      _kPlus_space[iState][iEq] = max(0.,_k[iState][iEq]);
      _kMin[iState][iEq]  = min(0.,_k[iState][iEq]*tCoeff + oneforthArea);
      _kPlus[iState][iEq] = max(0.,_k[iState][iEq]*tCoeff + _cellVolume/4.0);
      mini = min(mini, (Area+PastArea)/(4.0*_kPlus_space[iState][iEq]));

    }

    CFLogDebugMax( "kPlus = " << _kPlus[iState] << "\n");
    CFLogDebugMax( "kMin = " << _kMin[iState] << "\n");
    CFLogDebugMax( "_adimNormal = " << _adimNormal << "\n");
    CFLogDebugMax( "nodeArea = " << _nodeArea[iState] << "\n");

    if (_isOnlySplitter && !getMethodData().getDistributionData().isPerturb) {
      _updateCoeff[states[iState]->getLocalID()] += (_cellVolume+(tCoeff*_kPlus_space[iState].max()));

    }
  }

  if (!getMethodData().getDistributionData().isPerturb) {
    // Compute max timestep
    CFreal maxDT = SubSystemStatusStack::getActive()->getMaxDT();
    ///@todo this needs to be changed for hybrid meshes
    CFreal DT = mini/(DIM+1);
    if (DT < maxDT){
      SubSystemStatusStack::getActive()->setMaxDT(DT);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
