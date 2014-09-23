#include "MarcoTest/MarcoTest.hh"
#include "MarcoTest/StdUnSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/ProxyDofIterator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUnSetup,
                      MarcoTestMethodData,
                      MarcoTestModule>
stdUnSetupMarcoTestProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(const std::string& name) : 
  MarcoTestMethodCom(name),
  socket_gstates("gstates"),
  socket_nstatesProxy("nstatesProxy")
{
}

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::~StdUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstatesProxy);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;

  deleteAllPtr(socket_gstates);
  deleteAllPtr(socket_nstatesProxy);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MarcoTest

} // namespace COOLFluiD
