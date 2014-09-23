#include "Framework/MeshData.hh"
#include "Framework/JacobianLinearizer.hh"

#include "FluctSplit/ArtificialDiffusionStrategy.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

ArtificialDiffusionStrategy::ArtificialDiffusionStrategy(const std::string& name) :
  Framework::MethodStrategy<FluctuationSplitData>(name),
  m_result(0.,0),
  socket_normals("normals"),
  socket_isBState("isBState"),
  socket_volumes("volumes"),
  _solutionToDistMatTrans(CFNULL),
  _distToSolutionMatTrans(CFNULL),
  _linearToDistMatTrans(CFNULL),
  _solutionToLinearInUpdateMatTrans(CFNULL),
  _solutionToLinearMatTrans(CFNULL),
  _updateToLinearVecTrans(CFNULL),
  _linearizer(CFNULL),
  _distribVar(CFNULL),
  _phiS(),
  m_nodalST(),
  m_cell(CFNULL),
  m_states(CFNULL),
  m_tStates(CFNULL),
  m_linearStates(CFNULL),
  _isBoundaryState(),
  _stdTrsGeoBuilder(),
  m_cellID(0)
{
}
ArtificialDiffusionStrategy::~ArtificialDiffusionStrategy()
{
}
//////////////////////////////////////////////////////////////////////////////

void ArtificialDiffusionStrategy::configure ( Config::ConfigArgs& args )
{
  Framework::MethodStrategy<FluctuationSplitData>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ArtificialDiffusionStrategy::setup()
{
  CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
//   CFuint dim   = PhysicalModelStack::getActive()->getDim();

  _phiS.resize(nbEqs);

  m_nodalST.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());
  for (CFuint i =0; i < m_nodalST.size(); ++i) {
    m_nodalST[i].resize(nbEqs);
  }

  // get the linearizer
  _linearizer = getMethodData().getLinearizer();

  // get the distribution var set
  _distribVar = getMethodData().getDistribVar();

  _solutionToDistMatTrans = getMethodData().getSolutionToDistribMatTrans();
  _distToSolutionMatTrans = getMethodData().getDistribToSolutionMatTrans();
  _linearToDistMatTrans = getMethodData().getLinearToDistribMatTrans();
  _solutionToLinearInUpdateMatTrans =
    getMethodData().getSolutionToLinearInUpdateMatTrans();
  _solutionToLinearMatTrans = getMethodData().getSolutionToLinearMatTrans();
  _updateToLinearVecTrans = getMethodData().getUpdateToLinearVecTrans();

  // set up the geometric entity builder
  _stdTrsGeoBuilder.setup();

  _isBoundaryState.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
ArtificialDiffusionStrategy::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;

   result.push_back(&socket_normals);
   result.push_back(&socket_isBState);
   result.push_back(&socket_volumes);

   return result;
}
    } //namespace Fluctsplit

} //namespace COOLFluid
