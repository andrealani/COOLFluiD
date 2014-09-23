#include "RKRD/RungeKuttaRD.hh"
#include "StdUnSetup.hh"
#include "MathTools/RealVector.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUnSetup, RKRDData, RKRDModule> stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(const std::string& name) : RKRDCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_kstates("kstates"),
  socket_median_areas("median_areas")
{
}

void StdUnSetup::execute()
{
  CFAUTOTRACE;

  socket_rhs.getDataHandle().resize(0);

  socket_updateCoeff.getDataHandle().resize(0);

  socket_median_areas.getDataHandle().resize(0);

  socket_kstates.getDataHandle().resize(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdUnSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_median_areas);
  result.push_back(&socket_kstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD
