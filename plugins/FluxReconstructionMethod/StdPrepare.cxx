#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/StdPrepare.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdPrepare, FluxReconstructionSolverData, FluxReconstructionModule> stdPrepareProvider("StdPrepare");

//////////////////////////////////////////////////////////////////////////////

StdPrepare::StdPrepare(const std::string& name) : FluxReconstructionSolverCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  socket_gradientsAV("gradientsAV")
{
}

//////////////////////////////////////////////////////////////////////////////

StdPrepare::~StdPrepare()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdPrepare::execute()
{
  CFAUTOTRACE;

  // reset the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
  rhs = 0.0;

  // reset the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff = 0.0;

  // reset the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  const CFuint nGrads = gradients.size();
  for (CFuint iGrad = 0; iGrad < nGrads; ++iGrad)
  {
    const CFuint nVar = gradients[iGrad].size();
    for (CFuint iVar = 0; iVar < nVar; ++iVar)
    {
      gradients[iGrad][iVar] = 0.0;
    }
  }
  
  // reset the gradientsAV
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();
  const CFuint nGrads2 = gradientsAV.size();
  for (CFuint iGrad = 0; iGrad < nGrads2; ++iGrad)
  {
    const CFuint nVar = gradientsAV[iGrad].size();
    for (CFuint iVar = 0; iVar < nVar; ++iVar)
    {
      gradientsAV[iGrad][iVar] = 0.0;
    }
  }

  // reset the jacobian
  if (getMethodData().getLinearSystemSolver().size() > 0)
  {
    SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLinearSystemSolver()[0]->getMatrix();
    if (jacobMatrix.isNotNull())
    {
      jacobMatrix->resetToZeroEntries();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdPrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);
  result.push_back(&socket_gradientsAV);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
