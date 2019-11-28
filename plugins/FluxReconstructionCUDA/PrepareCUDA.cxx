#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionCUDA/FluxReconstructionCUDA.hh"
#include "FluxReconstructionCUDA/PrepareCUDA.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PrepareCUDA, FluxReconstructionSolverData, FluxReconstructionCUDAModule> prepareCUDAProvider("PrepareCUDA");

//////////////////////////////////////////////////////////////////////////////

PrepareCUDA::PrepareCUDA(const std::string& name) : StdPrepare(name),
  socket_gradientsCUDA("gradientsCUDA")
{
}

//////////////////////////////////////////////////////////////////////////////

PrepareCUDA::~PrepareCUDA()
{
}

//////////////////////////////////////////////////////////////////////////////

void PrepareCUDA::execute()
{
  CFAUTOTRACE;

  StdPrepare::execute();

  // reset the gradients
  DataHandle< CFreal > gradients = socket_gradientsCUDA.getDataHandle();
  const CFuint nGrads = gradients.size();
  for (CFuint iGrad = 0; iGrad < nGrads; ++iGrad)
  {
    gradients[iGrad] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
PrepareCUDA::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdPrepare::needsSockets();

  result.push_back(&socket_gradientsCUDA);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
