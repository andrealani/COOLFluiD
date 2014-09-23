#include "STM_SplitterSys.hh"
#include "Common/CFLog.hh"
#include "Framework/SubSystemStatus.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/MatrixInverter.hh"
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

STM_SplitterSys::STM_SplitterSys(const std::string& name) :
  SpaceTime_Splitter(name),
  _identity(),
  _eValues(0),
  _eValuesP(0),
  _eValuesM(0),
  _tempKp(),
  _tempKm(),
  _inverter(CFNULL),
  _kPlus(0),
  _kMin(0)
{
}

//////////////////////////////////////////////////////////////////////////////

STM_SplitterSys::~STM_SplitterSys()
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

void STM_SplitterSys::setBlockData()
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

void STM_SplitterSys::setup()
{
  SpaceTime_Splitter::setup();

  _identity.resize(PhysicalModelStack::getActive()->getNbEq());
  _identity =  0.00000001;

  CFLogDebugMax( "STM_SplitterSys::_nbEquations: "
  << _nbEquations
  << ", _firstVar = " << _firstVarID
  << ", _lastVar = " << _lastVarID
  << "\n");

  _eValues.resize(PhysicalModelStack::getActive()->getNbEq());
  _eValuesP.resize(PhysicalModelStack::getActive()->getNbEq());
  _eValuesM.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempKp.resize(_nbEquations,_nbEquations);
  _tempKm.resize(_nbEquations,_nbEquations);
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

}

//////////////////////////////////////////////////////////////////////////////

void STM_SplitterSys::computeK(const vector<State*>& states,
					const InwardNormalsData* const normalsData)
{
  m_normals = normalsData;
  _nbStatesInCell = states.size();
  _nodeArea.resize(states.size());


  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFreal kCoeff = 1./PhysicalModelStack::getActive()->getDim();
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1);
  CFreal mini = 10.e+10;
  CFreal temp;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

    CFLogDebugMax( "iState = " << iState << "\n");

    SpaceTime_Splitter::setAdimensionalNormal(iState);

    // if Moving Mesh, then substract the cell speed (on the diagonal)
    if (SubSystemStatusStack::getActive()->isMovingMesh() != false){

      // Here the tempKp and tempKm are used for storing
      // the right and left eigenvectors !!!!!!!!!
      getMethodData().getDistribVar()->computeEigenValuesVectors(_tempKp,
							     _tempKm,
							     _eValues,
							     _adimNormal);

      temp = 0.;
      //CFout << _cellSpeed << "\n";
      for (CFuint iDim = 0; iDim < nbDim; ++iDim){
	temp += _cellSpeed[iDim]*_adimNormal[iDim];
      }

      // Modify the Eigen Values to take into account the speed of the cell
      for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
	_eValues[iEq] -= temp;
	_eValuesP[iEq] = max(0.,_eValues[iEq]);
	_eValuesM[iEq] = min(0.,_eValues[iEq]);
      }

      // compute jacobian + and -
      *_kPlus[iState] = _tempKp*(_eValuesP*_tempKm);
      *_kMin[iState] = _tempKp*(_eValuesM*_tempKm);
    }
    else{
      getMethodData().getDistribVar()->splitJacobian(*_kPlus[iState],
						  *_kMin[iState],
						  _eValues,
						  _adimNormal);
    }

    _nodeArea[iState] = m_normals->getAreaNode(iState);

    *_kPlus[iState] *= kCoeff * _nodeArea[iState];
    *_kMin[iState]  *= kCoeff * _nodeArea[iState];
    _eValues *= kCoeff * _nodeArea[iState];
    // add Eps on the diagonal of K- and K+
    *_kMin[iState] -= _identity;
    *_kPlus[iState] += _identity;

    CFLogDebugMax( "kPlus @STM_SplitterSys::computeK" << "\n"
    << *_kPlus[iState] << "\n");
    CFLogDebugMax( "kMin  @STM_SplitterSys::computeK" << "\n"
    << *_kMin[iState]  << "\n");

    if (!getMethodData().getDistributionData().isPerturb) {
      // unused //    const CFreal maxEigenValue = max(0., _eValues.max());
      const CFreal maxKiplus = max(0., (*_kPlus[iState]).max());
      //CFout << maxKiplus << "  " << maxEigenValue << "\n";
      const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();

      mini = min(mini, (_cellVolume+_pastCellVolume)/(2*maxKiplus));

      _updateCoeff[states[iState]->getLocalID()] += (Area+(timeStep*maxKiplus*0.5));
      //_updateCoeff[states[iState]->getLocalID()] += Area/(Area+((timeStep/2.)*maxEigenValue));

      CFLogDebugMax( "updateCoeff @STM_SplitterSys::computeK" << "\n"
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
