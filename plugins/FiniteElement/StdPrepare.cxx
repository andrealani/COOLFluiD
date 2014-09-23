#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/StdPrepare.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdPrepare, FiniteElementMethodData, FiniteElementModule> stdPrepareProvider("StdPrepare");

//////////////////////////////////////////////////////////////////////////////

StdPrepare::StdPrepare(const std::string& name) : FiniteElementMethodCom(name),
  socket_isUpdated("isUpdated")
{
}

//////////////////////////////////////////////////////////////////////////////

void StdPrepare::execute()
{
  CFAUTOTRACE;

  if(!getMethodData().isSysMatrixFrozen()) {
    for (CFuint i = 0; i < getMethodData().getLinearSystemSolver().size(); ++i) {
      SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[i];
      SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

      // reset to zero all non zero entries in the jacobian
      jacobMatrix->resetToZeroEntries();
    }
  }

  // flag to check if a state has been updated or not
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  isUpdated = false;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdPrepare::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
