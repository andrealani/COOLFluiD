#include "STKS_SplitterSys.hh"
#include "Common/CFLog.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "MathTools/MatrixInverter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

STKS_SplitterSys::STKS_SplitterSys(const std::string& name) :
  SpaceTime_Splitter(name),
  _identity(),
  _eValues(0),
  _inverter(CFNULL),
  _kPlus(0),
  _kMin(0)
{
}

//////////////////////////////////////////////////////////////////////////////

STKS_SplitterSys::~STKS_SplitterSys()
{
  for (CFuint i = 0; i < _kPlus.size(); ++i) {
    deletePtr(_kPlus[i]);
  }

  for (CFuint i = 0; i < _kPlus.size(); ++i) {
    deletePtr(_kMin[i]);
  }

  deletePtr(_inverter);
}

//////////////////////////////////////////////////////////////////////////////

void STKS_SplitterSys::setBlockData()
{
  const CFuint totalNbEq = PhysicalModelStack::getActive()->getNbEq();

  CFLogDebugMax( "RDS_SplitterSys::setBlockData() => " <<
		 "blockSeparator = " << _blockSeparator << "\n" <<
  "totalNbEq = " << totalNbEq << "\n");

  cf_assert(_blockSeparator <= totalNbEq);

  _nbEquations = _blockSeparator;
  _firstVarID = 0;
  _lastVarID = _blockSeparator;

  cf_assert(_nbEquations == _lastVarID - _firstVarID);
  cf_assert(_lastVarID <= totalNbEq);
}

//////////////////////////////////////////////////////////////////////////////

void STKS_SplitterSys::setup()
{
  SpaceTime_Splitter::setup();

  _identity.resize(PhysicalModelStack::getActive()->getNbEq());
  _identity =  0.00000001;

  CFLogDebugMax( "SpaceTimeRDS_SplitterSys::_nbEquations: "
  << _nbEquations
  << ", _firstVar = " << _firstVarID
  << ", _lastVar = " << _lastVarID
  << "\n");

  _eValues.resize(PhysicalModelStack::getActive()->getNbEq());
  _inverter = MatrixInverter::create(_nbEquations, false);

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _kPlus.resize(maxNbStatesInCell);
  _kMin.resize(maxNbStatesInCell);

  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _kPlus[i] = new RealMatrix(_nbEquations,
                                        _nbEquations);
  }

  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _kMin[i] = new RealMatrix(_nbEquations,
                                       _nbEquations);
  }

   DIM = PhysicalModelStack::getActive()->getDim();
}

//////////////////////////////////////////////////////////////////////////////

void STKS_SplitterSys::computeK(const vector<State*>& states,
					   const InwardNormalsData* const normalsData)
{
 m_normals = normalsData;
  _nbStatesInCell = states.size();
  _nodeArea.resize(states.size());

  const CFreal kCoeff = 1./DIM;
  const CFreal Area = _cellVolume/(DIM+1.);
  CFreal mini = 10.e+10;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

    CFLogDebugMax( "iState = " << iState << "\n");

    SpaceTime_Splitter::setAdimensionalNormal(iState);


    getMethodData().getDistribVar()->splitJacobian(*_kPlus[iState],
						  *_kMin[iState],
						  _eValues,
						  _adimNormal);

    _nodeArea[iState] = m_normals->getAreaNode(iState);

    *_kPlus[iState] *= kCoeff * _nodeArea[iState];
    *_kMin[iState]  *= kCoeff * _nodeArea[iState];
    _eValues *= kCoeff * _nodeArea[iState];
    // add Eps on the diagonal of K- and K+
    *_kMin[iState] -= _identity;
    *_kPlus[iState] += _identity;

    CFLogDebugMax( "kPlus @SpaceTimeRDS_SplitterSys::computeK" << "\n"
    << *_kPlus[iState] << "\n");
    CFLogDebugMax( "kMin  @SpaceTimeRDS_SplitterSys::computeK" << "\n"
    << *_kMin[iState]  << "\n");

    if (!getMethodData().getDistributionData().isPerturb) {
      // unused //    const CFreal maxEigenValue = max(0., _eValues.max());
      const CFreal maxKiplus = std::max(0., (*_kPlus[iState]).max());
      //CFout << maxKiplus << "  " << maxEigenValue << "\n";
      const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();

      mini = min(mini, (_cellVolume+_pastCellVolume)/(2*maxKiplus));

      _updateCoeff[states[iState]->getLocalID()] += (Area+(timeStep*maxKiplus*0.5));
      //_updateCoeff[states[iState]->getLocalID()] += Area/(Area+((timeStep/2.)*maxEigenValue));

      CFLogDebugMax( "updateCoeff @SpaceTimeRDS_SplitterSys::computeK" << "\n"
		     << _updateCoeff[states[iState]->getLocalID()] << "\n" << "\n");
    }
  }

  if (!getMethodData().getDistributionData().isPerturb) {
    // Compute max timestep
    CFreal maxDT = SubSystemStatusStack::getActive()->getMaxDT();
    ///@todo this needs to be changed for hybrid meshes
    CFreal DT = mini/(PhysicalModelStack::getActive()->getDim()+1);
    if (DT < maxDT) {
      SubSystemStatusStack::getActive()->setMaxDT(DT);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
