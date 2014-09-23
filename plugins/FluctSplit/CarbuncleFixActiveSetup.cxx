#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/CarbuncleFixActiveSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CarbuncleFixActiveSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aCarbuncleFixActiveSetupProvider("CarbuncleFixActiveSetup");

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFixActiveSetup::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
}

//////////////////////////////////////////////////////////////////////////////

CarbuncleFixActiveSetup::CarbuncleFixActiveSetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_fix_active("fix_active")
{
  addConfigOptionsTo(this);
  m_maxsubelems = 1;
  setParameter("MaxNbSubElems",&m_maxsubelems);
}

//////////////////////////////////////////////////////////////////////////////

CarbuncleFixActiveSetup::~CarbuncleFixActiveSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
CarbuncleFixActiveSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_fix_active);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CarbuncleFixActiveSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFixActiveSetup::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFixActiveSetup::execute()
{
  CFAUTOTRACE;

  const CFuint nbcells =
    MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  DataHandle<CFreal> fix_active = socket_fix_active.getDataHandle();
  fix_active.resize( nbcells );
//  fix_active.resize( nbcells * m_maxsubelems );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

