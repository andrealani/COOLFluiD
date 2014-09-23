#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/ThetaSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ThetaSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aThetaSetupProvider("ThetaSetup");

//////////////////////////////////////////////////////////////////////////////

void ThetaSetup::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
}

//////////////////////////////////////////////////////////////////////////////

ThetaSetup::ThetaSetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_thetas("thetas")
{
  addConfigOptionsTo(this);
  m_maxsubelems = 1;
  setParameter("MaxNbSubElems",&m_maxsubelems);
}

//////////////////////////////////////////////////////////////////////////////

ThetaSetup::~ThetaSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ThetaSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_thetas);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ThetaSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ThetaSetup::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ThetaSetup::execute()
{
  CFAUTOTRACE;

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbcells =
    MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  DataHandle<CFreal> thetas = socket_thetas.getDataHandle();
  thetas.resize( nbcells * nbEqs * m_maxsubelems );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

