#include "MatrixStabilityMethodWriter/MatrixStabilityMethodWriter.hh"
#include "MatrixStabilityMethodWriter/StdUnSetup.hh"
#include "MathTools/RealVector.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<
                      StdUnSetup,
                      MatrixStabilityMethodData,
                      MatrixStabilityMethodWriterModule
                     >
stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(const std::string& name) :
  MatrixStabilityMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff")
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdUnSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs.resize(0);

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD
