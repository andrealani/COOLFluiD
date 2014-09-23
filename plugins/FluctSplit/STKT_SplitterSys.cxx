#include "STKT_SplitterSys.hh"
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

STKT_SplitterSys::STKT_SplitterSys(const std::string& name) :
  SpaceTime_Splitter(name),
  _identity(),
  _eValues(0),
  _eValuesP(0),
  _eValuesM(0),
  _rightEv(),
  _leftEv(),
  _inverter(CFNULL),
  _kPlus(0),
  _kMin(0)
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_SplitterSys::~STKT_SplitterSys()
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

void STKT_SplitterSys::setBlockData()
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

void STKT_SplitterSys::setup()
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
  _eValuesP.resize(PhysicalModelStack::getActive()->getNbEq());
  _eValuesM.resize(PhysicalModelStack::getActive()->getNbEq());
  _rightEv.resize(_nbEquations,_nbEquations);
  _leftEv.resize(_nbEquations,_nbEquations);
  _inverter = MatrixInverter::create(_nbEquations, false);

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _kPlus.resize(maxNbStatesInCell);
  _kMin.resize(maxNbStatesInCell);
  // We assume that we have only one type of element
  // Then maxNbStatesInCell is also the number of states per cell
_nbStatesInCell=maxNbStatesInCell;
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

void STKT_SplitterSys::computeK(const vector<State*>& states,
					   const InwardNormalsData* const normalsData)
{
  m_normals = normalsData;
  _nbStatesInCell = states.size();
  
  const CFreal kCoeff = 1./DIM;
  CFreal Area = _cellVolume;
  CFreal PastArea = _pastCellVolume;
  const CFreal tCoeff = _timeStep/2.;
  CFreal mini = 10.e+10;
  CFreal onetheirdArea = _cellVolume/(DIM+1.0);
  CFreal knodeArea;
  DistributionData& ddata = getMethodData().getDistributionData();
  ///@TODO Since for the moment the High order is not working, I have commented
  /// this part !!!!!!!!
  //THe following is valid only fot P2 in 2D
   if (getMethodData().getDistributionData().isHO){
       Area /= 4.0;
       PastArea /= 4.0;
       onetheirdArea /= 4.0;
   }
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

    knodeArea = (m_normals->getAreaNode(iState))*kCoeff;
    CFLogDebugMax( "iState = " << iState << "\n");

    SpaceTime_Splitter::setAdimensionalNormal(iState);

    getMethodData().getDistribVar()->computeEigenValuesVectors(_rightEv, _leftEv, _eValues,_adimNormal);

    // First we use the positive part of eigen value to compute dt and the update coeff
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {

	_eValuesP[iEq] = knodeArea*max(0.,_eValues[iEq]);
      }

    if (!ddata.isPerturb) {

      const CFreal maxKiplus = max(0.0, _eValuesP.max());
      if (maxKiplus != 0.0)
        mini = min(mini, (Area+PastArea)/(maxKiplus));
      _updateCoeff[states[iState]->getLocalID()] += (onetheirdArea+(tCoeff*maxKiplus));
      CFLogDebugMax( "updateCoeff @SpaceTimeRDS_SplitterSys::computeK" << "\n"
		     << _updateCoeff[states[iState]->getLocalID()] << "\n" << "\n");
    }
    
   //Modification of the eigen value to comupte the space time k
   knodeArea *= tCoeff;
   for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
	_eValues[iEq] = knodeArea*_eValues[iEq] + onetheirdArea;
	_eValuesP[iEq] = max(0.,_eValues[iEq]);
	_eValuesM[iEq] = min(0.,_eValues[iEq]);
      }

    // compute jacobian + and -
      *_kPlus[iState] = _rightEv*(_eValuesP*_leftEv);
      *_kMin[iState] = _rightEv*(_eValuesM*_leftEv);

   CFLogDebugMax( "kPlus @SpaceTimeRDS_SplitterSys::computeK" << "\n"
    << *_kPlus[iState] << "\n");
    CFLogDebugMax( "kMin  @SpaceTimeRDS_SplitterSys::computeK" << "\n"
    << *_kMin[iState]  << "\n");
  }

  if (!getMethodData().getDistributionData().isPerturb) {
    // Compute max timestep
    CFreal maxDT = SubSystemStatusStack::getActive()->getMaxDT();
    ///@todo this needs to be changed for hybrid meshes
    ///@todo this needs to be changed for hybrid meshes
    CFreal DT = mini/(DIM+1);
    if (DT < maxDT) {
      SubSystemStatusStack::getActive()->setMaxDT(DT);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
