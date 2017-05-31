#include "FluctSplit/Null_ArtDiffStrategy.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Null_ArtDiffStrategy,
                       FluctuationSplitData,
                       ArtificialDiffusionStrategy,
                       FluctSplitModule>
                       nullArtDiffStrategyProvider("Null");

//////////////////////////////////////////////////////////////////////////////

Null_ArtDiffStrategy::Null_ArtDiffStrategy(const std::string& name) :
  ArtificialDiffusionStrategy(name)
{
}

//////////////////////////////////////////////////////////////////////////////

Null_ArtDiffStrategy::~Null_ArtDiffStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void Null_ArtDiffStrategy::setup()
{
  ArtificialDiffusionStrategy::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& Null_ArtDiffStrategy::operator()(const std::vector<Framework::State*>& vars,
                            const RealVector& shapeF,
                            const RealMatrix& grad,
                            Framework::GeometricEntity* const geo)
{
  // this is just not to have a warning
  static RealVector avector;
  return avector;
}

//////////////////////////////////////////////////////////////////////////////

void Null_ArtDiffStrategy::addArtificialDiff(std::vector<RealVector>& residual)
{
}

//////////////////////////////////////////////////////////////////////////////

    }// End namespace FluctSplit

}// End namespace COOLFluiD
