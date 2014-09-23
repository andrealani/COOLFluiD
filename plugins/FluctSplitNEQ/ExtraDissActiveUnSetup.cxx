#include "Framework/MethodCommandProvider.hh"

#include "ExtraDissActiveUnSetup.hh"
#include "FluctSplitNEQ.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ExtraDissActiveUnSetup,
                      FluctSplit::FluctuationSplitData,
                      FluctSplitNEQModule>
anExtraDissActiveUnSetupProvider("ExtraDissActiveUnSetup");

//////////////////////////////////////////////////////////////////////////////

ExtraDissActiveUnSetup::ExtraDissActiveUnSetup(const std::string& name) :
  FluctSplit::FluctuationSplitCom(name),
  socket_ExtraDiss_active("extra_dissipation_active")
{
}

//////////////////////////////////////////////////////////////////////////////

ExtraDissActiveUnSetup::~ExtraDissActiveUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ExtraDissActiveUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_ExtraDiss_active);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ExtraDissActiveUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal> ExtraDiss_active = socket_ExtraDiss_active.getDataHandle();
  ExtraDiss_active.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ

} // namespace COOLFluiD
