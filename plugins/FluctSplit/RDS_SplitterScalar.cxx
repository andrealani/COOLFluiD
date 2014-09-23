#include "FluctSplit/RDS_SplitterScalar.hh"
#include "Common/CFLog.hh"
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

RDS_SplitterScalar::RDS_SplitterScalar(const std::string& name) :
  Splitter(name),
  socket_updateCoeff("updateCoeff"),
  _kPlus(0),
  _kMin(0),
  _k(0)
{
}

//////////////////////////////////////////////////////////////////////////////

RDS_SplitterScalar::~RDS_SplitterScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterScalar::setBlockData()
{
  const CFuint totalNbEq = PhysicalModelStack::getActive()->getNbEq();

  cf_assert(_blockSeparator <= totalNbEq);

  if (_blockSeparator == totalNbEq)
  {
    _firstVarID = 0;
    _lastVarID  = totalNbEq;
    _nbEquations = totalNbEq;
  }
  else
  {
    if (!getMethodData().isScalarFirst()) {
      _firstVarID = _blockSeparator;
      _lastVarID = totalNbEq;
      _nbEquations = _lastVarID - _firstVarID;
    }
    else {
      _nbEquations = _blockSeparator;
      _firstVarID = 0;
      _lastVarID = _blockSeparator; 
    }
  }
  
  cf_assert(_nbEquations == _lastVarID - _firstVarID);
  cf_assert(_lastVarID <= totalNbEq);
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterScalar::setup()
{
  Splitter::setup();

  _isOnlySplitter =
     (_nbEquations == PhysicalModelStack::getActive()->getNbEq());

  m_invDim = 1.0 / PhysicalModelStack::getActive()->getDim();

  m_maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  _kPlus.resize(m_maxNbStatesInCell);
  _kMin.resize(m_maxNbStatesInCell);
  _k.resize(m_maxNbStatesInCell);

  for (CFuint i = 0; i < m_maxNbStatesInCell; ++i) {
    _kPlus[i].resize(_nbEquations);
    _kMin[i].resize(_nbEquations);
    _k[i].resize(_nbEquations);
  }
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterScalar::doComputeK(CFuint iState)
{
  // The transformation of the normal is needed if the coordinate system
  // is rotated. The normal is adimensionalized, so there is need to multiply
  // by the nodeArea when computing the k parameter

  Splitter::setAdimensionalNormal(iState);

  // uses the linearized state to compute the K vector
  getMethodData().getDistribVar()->computeScalarJacobian(_adimNormal, _k[iState]);

  CFreal nodeArea = m_normals->getAreaNode(iState);
  nodeArea   *= m_invDim;
  _k[iState] *= nodeArea;

  // calculate Kplus and Kminus
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq)
  {
    _kPlus[iState][iEq] = std::max(0.,_k[iState][iEq]);
    _kMin[iState][iEq]  = std::min(0.,_k[iState][iEq]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterScalar::computeK(const vector<State*>& states,
				  const InwardNormalsData* const normalsData)
{
  m_normals = normalsData;
  _nbStatesInCell = states.size();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // The transformation of the normal is needed if the coordinate system is rotated
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    doComputeK(iState);

    if (_isOnlySplitter && !getMethodData().getDistributionData().isPerturb) {
      updateCoeff[states[iState]->getLocalID()] += _kPlus[iState].max();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
RDS_SplitterScalar::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
    Splitter::needsSockets();

  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

