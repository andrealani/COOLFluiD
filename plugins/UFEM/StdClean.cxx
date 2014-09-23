#include "Framework/MethodCommandProvider.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/StdClean.hh"
#include "Framework/LSSMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdClean, UFEMSolverData, UFEMPlugin> StdCleanProvider("StdClean");

//////////////////////////////////////////////////////////////////////////////

StdClean::StdClean(const std::string& name) :
  UFEMSolverCom(name),
  socket_isUpdated("isUpdated"),
  socket_rhs("rhs")
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

StdClean::~StdClean()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void StdClean::execute()
{
  CFAUTOTRACE;

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated = false;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs = 0.;

  SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLinearSystemSolver()[0]->getMatrix();
  jacobMatrix->resetToZeroEntries();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdClean::needsSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_rhs);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace UFEM

} // namespace COOLFluiD
