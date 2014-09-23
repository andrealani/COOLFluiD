#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/CarbuncleFixActiveUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CarbuncleFixActiveUnSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aCarbuncleFixActiveUnSetupProvider("CarbuncleFixActiveUnSetup");

//////////////////////////////////////////////////////////////////////////////

CarbuncleFixActiveUnSetup::CarbuncleFixActiveUnSetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_fix_active("fix_active")
{
}

//////////////////////////////////////////////////////////////////////////////

CarbuncleFixActiveUnSetup::~CarbuncleFixActiveUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CarbuncleFixActiveUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_fix_active);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFixActiveUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal> fix_active = socket_fix_active.getDataHandle();
  fix_active.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
