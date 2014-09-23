#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/UpdateCoeffSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateCoeffSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aUpdateCoeffSetupProvider("UpdateCoeffSetup");

//////////////////////////////////////////////////////////////////////////////

void UpdateCoeffSetup::defineConfigOptions(Config::OptionList& options)
{
//    options.addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
}

//////////////////////////////////////////////////////////////////////////////

UpdateCoeffSetup::UpdateCoeffSetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_updateCoeff("updateCoeff")
{
  addConfigOptionsTo(this);
//   m_maxsubelems = 1;
//   setParameter("MaxNbSubElems",&m_maxsubelems);
}

//////////////////////////////////////////////////////////////////////////////

UpdateCoeffSetup::~UpdateCoeffSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
UpdateCoeffSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpdateCoeffSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UpdateCoeffSetup::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void UpdateCoeffSetup::execute()
{
  CFAUTOTRACE;

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbcells =
    MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize( nbcells );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

