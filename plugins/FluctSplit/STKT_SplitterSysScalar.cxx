#include "STKT_SplitterSysScalar.hh"
#include "Common/CFLog.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "NavierStokes/Euler2DSymm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////
      
void STKT_SplitterSysScalar::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< CFuint >("SysSize","Size of the coupled system of equations.");
  options.addConfigOption< CFuint >("SysStartID","ID corresponding to the start of the system.");
  options.addConfigOption< vector<CFuint> >("ScalarEqIDs","IDs corresponding to the decoupled scalar equations.");
}
      
//////////////////////////////////////////////////////////////////////////////

STKT_SplitterSysScalar::STKT_SplitterSysScalar(const std::string& name) :
  STKT_SplitterSys(name),
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
      
STKT_SplitterSysScalar::~STKT_SplitterSysScalar()
{
}

//////////////////////////////////////////////////////////////////////////////
      
void STKT_SplitterSysScalar::setBlockData()
{
  CFLogDebugMin( "STKT_SplitterSys::setBlockData() => " <<
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
      
void STKT_SplitterSysScalar::setup()
{
  CFAUTOTRACE;
  
  STKT_SplitterSys::setup(); 
  
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _kPlusScalar.resize(maxNbStatesInCell);
  _kMinScalar.resize(maxNbStatesInCell);
  _kScalar.resize(maxNbStatesInCell);
  
  const CFuint scalarSize = _scalarEqIDs.size();
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _kPlusScalar[i].resize(scalarSize);
    _kMinScalar[i].resize(scalarSize);
    _kScalar[i].resize(scalarSize);
  }
}

//////////////////////////////////////////////////////////////////////////////

void STKT_SplitterSysScalar::computeK(const std::vector<Framework::State*>& states,
				      const InwardNormalsData* const normalsData)
{
  m_normals = normalsData;
  _nbStatesInCell = states.size();
  
  const CFreal kCoeff = 1./DIM;
  CFreal Area = _cellVolume;
  CFreal PastArea = _pastCellVolume;
  const CFreal tCoeff = _timeStep/2.;
  CFreal mini = 10.e+10;
  _oneThirdVolume = _cellVolume/(DIM+1.0);
  _knodeVolume = 0.;
  DistributionData& ddata = getMethodData().getDistributionData();
  
  using namespace Physics::NavierStokes;
  static SafePtr<Euler2DSymm> dvar = getMethodData().getDistribVar().d_castTo<Euler2DSymm>();
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _knodeVolume = (m_normals->getAreaNode(iState))*kCoeff*tCoeff;
    
    Splitter::setAdimensionalNormal(iState);
    
    dvar->splitJacobian(*_kPlus[iState],*_kMin[iState], _eValues,_adimNormal, *this);
    
    // First we use the positive part of eigen value to compute dt and the update coeff
    for (CFuint iEq = 0; iEq < _sysSize; ++iEq) {
      _eValuesP[iEq] = _knodeVolume*max(0.,_eValues[iEq]);
    }
    
    if (!ddata.isPerturb) {
      const CFreal maxKiplus = max(0.0, _eValuesPmax());
      if (maxKiplus > 0.0) {
        mini = min(mini, (Area+PastArea)/maxKiplus);
      }
      _updateCoeff[states[iState]->getLocalID()] += (_oneThirdVolume + maxKiplus);
      CFLogDebugMax( "updateCoeff @SpaceTimeRDS_SplitterSys::computeK" << "\n"
		     << _updateCoeff[states[iState]->getLocalID()] << "\n" << "\n");
    }
    
    // uses the linearized state to compute the K vector
    dvar->computeScalarJacobian(_adimNormal, _kScalar[iState]);
    
    // calculate Kplus and Kminus
    const CFuint scalarSize = _scalarEqIDs.size();
    for (CFuint iEq = 0; iEq < scalarSize; ++iEq) { 
      _kScalar[iState][iEq] = _oneThirdVolume + _knodeVolume*_kScalar[iState][iEq];
      _kPlusScalar[iState][iEq] = max(0.,_kScalar[iState][iEq]);
      _kMinScalar[iState][iEq]  = min(0.,_kScalar[iState][iEq]);
    }
    
//     RealSliceMatrix::setNbRowsCols(3,3);
//     RealMatrix kp(4,4);
//     RealMatrix km(4,4);
//     kp.slice(0,0) = *_kPlus[iState];  kp(3,3) = _kPlusScalar[iState][0];
//     km.slice(0,0) = *_kMin[iState];   km(3,3) = _kMinScalar[iState][0];
//     cout << "kp[" << iState <<"]" << endl; cout.precision(14); cout << kp << endl;
//     cout << "km[" << iState <<"]" << endl; cout.precision(14); cout << km << endl;
  } 
  
 
  

  
  if (!ddata.isPerturb) {
    // Compute max timestep
    const CFreal maxDT = SubSystemStatusStack::getActive()->getMaxDT();
    ///@todo this needs to be changed for hybrid meshes
    const CFreal DT = mini/(DIM+1);
    if (DT < maxDT) {
      SubSystemStatusStack::getActive()->setMaxDT(DT);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
