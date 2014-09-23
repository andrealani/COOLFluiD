#include "RDS_SplitterSysScalar.hh"
#include "Common/CFLog.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////
      
void RDS_SplitterSysScalar::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< CFuint >("SysSize","Size of the coupled system of equations.");
  options.addConfigOption< CFuint >("SysStartID","ID corresponding to the start of the system.");
  options.addConfigOption< vector<CFuint> >("ScalarEqIDs","IDs corresponding to the decoupled scalar equations.");
}
      
//////////////////////////////////////////////////////////////////////////////

RDS_SplitterSysScalar::RDS_SplitterSysScalar(const std::string& name) :
  RDS_SplitterSys(name),
  _kPlusScalar(0),
  _kMinScalar(0),
  _kScalar(0)
{ 
  this->addConfigOptionsTo(this);
  
  _sysSize = 0;
  this->setParameter("SysSize",&_sysSize);
  
  _sysStartID = 0;
  this->setParameter("SysStartID",&_sysStartID);
  
  _scalarEqIDs = vector<CFuint>();
  this->setParameter("ScalarEqIDs",&_scalarEqIDs);
}
      
//////////////////////////////////////////////////////////////////////////////
      
RDS_SplitterSysScalar::~RDS_SplitterSysScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSysScalar::setBlockData()
{
  CFLogDebugMin( "RDS_SplitterSys::setBlockData() => " <<
		 "blockSeparator = " << _blockSeparator << "\n" <<
		 "totalNbEq = " <<  PhysicalModelStack::getActive()->getNbEq() << "\n");
  
  cf_assert(_blockSeparator <=  PhysicalModelStack::getActive()->getNbEq() );
  
  // special treatement: this applies to system only
  _nbEquations = _sysSize;
  _firstVarID = _sysStartID;
  _lastVarID = _sysStartID + _sysSize;
  
  cf_assert(_nbEquations == _lastVarID - _firstVarID);
  cf_assert(_lastVarID <= PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////
      
void RDS_SplitterSysScalar::setup()
{
  CFAUTOTRACE;
  
  RDS_SplitterSys::setup(); 
  
  _kPlusScalar.resize(m_maxNbStatesInCell);
  _kMinScalar.resize(m_maxNbStatesInCell);
  _kScalar.resize(m_maxNbStatesInCell);
  
  const CFuint scalarSize = _scalarEqIDs.size();
  if (scalarSize > 0) {
    for (CFuint i = 0; i < m_maxNbStatesInCell; ++i) {
      _kPlusScalar[i].resize(scalarSize);
      _kMinScalar[i].resize(scalarSize);
      _kScalar[i].resize(scalarSize);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSysScalar::computeK(const vector<State*>& states,
				     const InwardNormalsData* const normalsData)
{
  m_normals = normalsData;
  _nbStatesInCell = states.size();
  
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  // apply the entropy or carbuncle fix
  getMethodData().getJacobianFixComputer()->computeFix(*normalsData, m_delta);
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    getMethodData().getDistribVar()->setDelta(m_delta[iState]);
    doComputeK(iState);
    
    if (!getMethodData().getDistributionData().isPerturb) {
      if (_sysSize > 0) {
	const CFreal maxEigenValue = std::max(0.0, m_eValues[iState]->max());
	updateCoeff[states[iState]->getLocalID()] += m_kCoeff*m_nodeArea*maxEigenValue;
      }
      else {
	updateCoeff[states[iState]->getLocalID()] += _kPlusScalar[iState].max();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSysScalar::doComputeK(CFuint iState)
{
  Splitter::setAdimensionalNormal(iState);
  
  m_nodeArea = m_normals->getAreaNode(iState);
  
  if (_sysSize > 0) {
    getMethodData().getDistribVar()->splitJacobian(*_kPlus[iState],
						   *_kMin[iState],
						   *m_eValues[iState],
						   _adimNormal);
    // system settings
    *_kPlus[iState] *= m_kCoeff * m_nodeArea;
    *_kMin[iState]  *= m_kCoeff * m_nodeArea;
  }
  
  const CFuint scalarSize = _scalarEqIDs.size();
  if (scalarSize > 0) { 
    // uses the linearized state to compute the K vector
    getMethodData().getDistribVar()->computeScalarJacobian(_adimNormal, _kScalar[iState]);
    
    // scalar settings
    _kScalar[iState] *= m_kCoeff * m_nodeArea;
    
    // calculate Kplus and Kminus
    for (CFuint iEq = 0; iEq < scalarSize; ++iEq) {
      _kPlusScalar[iState][iEq] = max(0.,_kScalar[iState][iEq]);
      _kMinScalar[iState][iEq]  = min(0.,_kScalar[iState][iEq]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
