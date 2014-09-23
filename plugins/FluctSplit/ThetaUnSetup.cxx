#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/ThetaUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ThetaUnSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aThetaUnSetupProvider("ThetaUnSetup");

//////////////////////////////////////////////////////////////////////////////

ThetaUnSetup::ThetaUnSetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_thetas("thetas")
{
}

//////////////////////////////////////////////////////////////////////////////

ThetaUnSetup::~ThetaUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ThetaUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_thetas);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ThetaUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal> thetas = socket_thetas.getDataHandle();
  thetas.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
