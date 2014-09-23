#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_StdComputeTimeRhsCoupling.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSIdxMapping.hh"
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

MethodCommandProvider<FVMCC_StdComputeTimeRhsCoupling, CellCenterFVMData, FiniteVolumeModule> fvmccStdComputeTimeRhsCoupling("StdTimeRhsCoupling");

//////////////////////////////////////////////////////////////////////////////

void FVMCC_StdComputeTimeRhsCoupling::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("useGlobalDT","Flag telling if to use global DT.");
  options.addConfigOption< bool >("useAnalyticalMatrix","Flag telling if to use analytical matrix.");
  options.addConfigOption< vector<bool> >
    ("annullDiagValue","Annull the diagonal value in a subsystem time contribution to the matrix.");
  options.addConfigOption< vector<bool> >
    ("zeroDiagValue", "Flags telling indicating entries in the matrix diagonal.");
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_StdComputeTimeRhsCoupling::FVMCC_StdComputeTimeRhsCoupling(const std::string& name) :
  CellCenterFVMCom(name),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  _lss(),
  _equations(),
  _jacobMatrix(),
  _idxMapping()
{
  addConfigOptionsTo(this);
  
  _useGlobalDT = false;
  setParameter("useGlobalDT",&_useGlobalDT);
  
  _useAnalyticalMatrix = false;
  setParameter("useAnalyticalMatrix",&_useAnalyticalMatrix);
  
  _annullDiagValue = vector<bool>();
  setParameter("annullDiagValue", &_annullDiagValue);

  _zeroDiagValue = vector<bool>();
  setParameter("zeroDiagValue",&_zeroDiagValue);
}
      
//////////////////////////////////////////////////////////////////////////////

FVMCC_StdComputeTimeRhsCoupling::~FVMCC_StdComputeTimeRhsCoupling()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_StdComputeTimeRhsCoupling::setup()
{
  CellCenterFVMCom::setup();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (_zeroDiagValue.size() != nbEqs) {
    CFLog(VERBOSE, "FVMCC_StdComputeTimeRhsCoupling::setup() => default values for _zeroDiagValue\n");
    _zeroDiagValue.resize(nbEqs);
    _zeroDiagValue.assign(nbEqs,false);
  }
  
  const CFuint nbLSS = getMethodData().getLinearSystemSolver().size();

  // linear system solver
  // acquaintance of the linear systems
  _lss.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    _lss[i] = getMethodData().getLinearSystemSolver()[i];
  }

  _equations.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    // acquaintance of the equation IDs to solve in each LSS
    _equations[i] = _lss[i]->getEquationIDs();
    cf_assert(_equations[i]->size() == _lss[i]->getNbSysEqs());
  }

  // acquaintance of the jacobian matrices
  _jacobMatrix.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    _jacobMatrix[i] = _lss[i]->getMatrix();
  }

  // acquaintance of the index mappings from local to global
  _idxMapping.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    _idxMapping[i] = &_lss[i]->getLocalToGlobalMapping();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_StdComputeTimeRhsCoupling::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbStates = states.size();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFuint nbLSS = _lss.size();

  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    if (states[iState]->isParUpdatable()) {

      // loop over the LSSs
      for (CFuint iLSS = 0; iLSS < nbLSS; ++iLSS) {
	const CFuint nbSysEqs = _equations[iLSS]->size();
	const CFreal diagValue = updateCoeff[iState*nbLSS + iLSS]/cfl;
	CFuint globalID = _idxMapping[iLSS]->getColID
	  (states[iState]->getLocalID())*nbSysEqs;

	for (CFuint iEq = 0; iEq < nbSysEqs; ++iEq, ++globalID) {
	  _jacobMatrix[iLSS]->addValue(globalID, globalID, diagValue);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_StdComputeTimeRhsCoupling::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
