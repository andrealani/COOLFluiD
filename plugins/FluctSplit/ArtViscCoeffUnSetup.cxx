#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/ArtViscCoeffUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ArtViscCoeffUnSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aArtViscCoeffUnSetupProvider("ArtViscCoeffUnSetup");

//////////////////////////////////////////////////////////////////////////////

ArtViscCoeffUnSetup::ArtViscCoeffUnSetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_ArtViscCoeff("ArtViscCoeff")
{
}

//////////////////////////////////////////////////////////////////////////////

ArtViscCoeffUnSetup::~ArtViscCoeffUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ArtViscCoeffUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_ArtViscCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ArtViscCoeffUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal> ArtViscCoeff = socket_ArtViscCoeff.getDataHandle();
  ArtViscCoeff.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
