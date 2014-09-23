#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_PseudoSteadyTimeRhsSingleSys.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_PseudoSteadyTimeRhsSingleSys,
		      CellCenterFVMData,
		      FiniteVolumeModule>
fvmccPseudoSteadyTimeRhsSingleSys("PseudoSteadyTimeRhsSingleSys");

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsSingleSys::FVMCC_PseudoSteadyTimeRhsSingleSys
(const std::string& name) :
  FVMCC_PseudoSteadyTimeRhsCoupling(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_PseudoSteadyTimeRhsSingleSys::~FVMCC_PseudoSteadyTimeRhsSingleSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PseudoSteadyTimeRhsSingleSys::execute()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbStates = states.size();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  // set the equations subsystem description
  const EquationSubSysDescriptor& eqSS =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint iLSS = eqSS.getEqSS();
  const vector<CFuint>& currEqs = *_equations[iLSS];
  _nbEqs = currEqs.size();
  _start = currEqs[0];
  const CFuint nbLSS = _lss.size();

  CFreal minDt = 0.0;
  if (_useGlobalDT) {
    // select the minimum delta T
    minDt = volumes[0]/updateCoeff[0];
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      minDt = min(volumes[iState]/updateCoeff[iState*nbLSS + iLSS], minDt);
    }
  }

  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {

    if (states[iState]->isParUpdatable()) {
      if (!_useGlobalDT) {
	_diagValue = (dt > 0.0) ? volumes[iState]/dt :
	  updateCoeff[iState*nbLSS + iLSS]/cfl;
      }
      else {
	_diagValue = volumes[iState]/(minDt*cfl);
      }

      if (!_useAnalyticalMatrix) {
	// compute the transformation matrix numerically
	computeNumericalTransMatrix(iState, iLSS);
      }
      else {
	computeAnalyticalTransMatrix(iState, iLSS);
      }

      // add the values in the jacobian matrix
      _jacobMatrix[iLSS]->addValues(*_acc[iLSS]);

      // reset to zero the entries in the block accumulator
      _acc[iLSS]->reset();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
