#include "FiniteVolume/FiniteVolume.hh"
#include "FVMCC_StdComputeTimeRhs.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

#ifdef CF_HAVE_CUDA
#include "Framework/CudaTimer.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_StdComputeTimeRhs, CellCenterFVMData, FiniteVolumeModule> fvmccStdComputeTimeRhs("StdTimeRhs");

//////////////////////////////////////////////////////////////////////////////

FVMCC_StdComputeTimeRhs::FVMCC_StdComputeTimeRhs(const std::string& name) :
  CellCenterFVMCom(name),
  _lss(CFNULL),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_rhs("rhs") 
{
  addConfigOptionsTo(this);
  
  _zeroDiagValue = vector<bool>();
  setParameter("zeroDiagValue",&_zeroDiagValue);
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_StdComputeTimeRhs::~FVMCC_StdComputeTimeRhs()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_StdComputeTimeRhs::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<bool> >
    ("zeroDiagValue", "Flags telling indicating entries in the matrix diagonal.");
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_StdComputeTimeRhs::setup()
{
  CellCenterFVMCom::setup();

  // linear system solver
  _lss = getMethodData().getLinearSystemSolver()[0];
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (_zeroDiagValue.size() != nbEqs) {
    CFLog(VERBOSE, "FVMCC_StdComputeTimeRhs::setup() => defualt values for _zeroDiagValue\n");
    _zeroDiagValue.resize(nbEqs);
    _zeroDiagValue.assign(nbEqs,false);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_StdComputeTimeRhs::execute()
{
  CFAUTOTRACE;

#ifdef CF_HAVE_CUDA
  CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
  timer.start();
#endif

  SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLSSMatrix(0);
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();

  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    if (states[iState]->isParUpdatable()) {

      const CFreal diagValue = updateCoeff[iState]/cfl;
      //CFuint globalID = states[iState]->getGlobalID()*nbEqs;
      CFuint globalID = idxMapping.getColID
	(states[iState]->getLocalID())*nbEqs;

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq, ++globalID) {
	jacobMatrix->addValue(globalID, globalID, diagValue);
      }
    }
  }
  
#ifdef CF_HAVE_CUDA  
  CFLog(VERBOSE, "FVMCC_StdComputeTimeRhs::execute() took " << timer.elapsed() << " s\n");
#endif
}
      
//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_StdComputeTimeRhs::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_rhs);
 
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
