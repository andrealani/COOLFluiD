#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/ArtViscCoeffSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ArtViscCoeffSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aArtViscCoeffSetupProvider("ArtViscCoeffSetup");

//////////////////////////////////////////////////////////////////////////////

void ArtViscCoeffSetup::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
}

//////////////////////////////////////////////////////////////////////////////

ArtViscCoeffSetup::ArtViscCoeffSetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_ArtViscCoeff("ArtViscCoeff")
{
  addConfigOptionsTo(this);
  m_maxsubelems = 1;
  setParameter("MaxNbSubElems",&m_maxsubelems);
}

//////////////////////////////////////////////////////////////////////////////

ArtViscCoeffSetup::~ArtViscCoeffSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ArtViscCoeffSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_ArtViscCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ArtViscCoeffSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ArtViscCoeffSetup::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ArtViscCoeffSetup::execute()
{
  CFAUTOTRACE;

  const CFuint nbcells =
   MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  DataHandle<CFreal> ArtViscCoeff = socket_ArtViscCoeff.getDataHandle();
  ArtViscCoeff.resize( nbcells );
  //  ArtViscCoeff.resize( nbcells * m_maxsubelems );????
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

