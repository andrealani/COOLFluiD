#include "Framework/MethodCommandProvider.hh"

#include "ExtraDissActiveSetup.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplitNEQ.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ExtraDissActiveSetup,
                      FluctSplit::FluctuationSplitData,
                      FluctSplitNEQModule>
anExtraDissActiveSetupProvider("ExtraDissActiveSetup");

//////////////////////////////////////////////////////////////////////////////

ExtraDissActiveSetup::ExtraDissActiveSetup(const std::string& name) :
  FluctSplit::FluctuationSplitCom(name),
  socket_ExtraDiss_active("extra_dissipation_active")
{
}

//////////////////////////////////////////////////////////////////////////////

ExtraDissActiveSetup::~ExtraDissActiveSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ExtraDissActiveSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_ExtraDiss_active);
  return result;
}


//////////////////////////////////////////////////////////////////////////////

void ExtraDissActiveSetup::configure ( Config::ConfigArgs& args )
{
  FluctSplit::FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ExtraDissActiveSetup::execute()
{
  CFAUTOTRACE;

  const CFuint nbcells =
    MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  DataHandle<CFreal> extraDiss_active = socket_ExtraDiss_active.getDataHandle();
  extraDiss_active.resize( nbcells );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ

} // namespace COOLFluiD

