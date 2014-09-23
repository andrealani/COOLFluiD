#include "Framework/MeshData.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/MeshMovementDirichletBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MeshMovementDirichletBC, FiniteElementMethodData, FiniteElementModule> MeshMovementDirichletBCProvider("MeshMovementDirichletBC");

//////////////////////////////////////////////////////////////////////////////

void MeshMovementDirichletBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("Implicit","Apply the BC implicitly?");
}

//////////////////////////////////////////////////////////////////////////////

MeshMovementDirichletBC::MeshMovementDirichletBC(const std::string& name) :
  FiniteElementMethodCom(name),
  _sockets(),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_states("states"),
  socket_bStatesNeighbors("bStatesNeighbors")
{
   addConfigOptionsTo(this);

  _isImplicit = false;
   setParameter("Implicit",&_isImplicit);

}

//////////////////////////////////////////////////////////////////////////////

MeshMovementDirichletBC::~MeshMovementDirichletBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshMovementDirichletBC::setup()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshMovementDirichletBC::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshMovementDirichletBC::executeOnTrs()
{
  CFAUTOTRACE;

  getMethodData().setDirichletBCApplied(true);

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "MeshMovementDirichletBC::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<std::valarray<Framework::State*> > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();

  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // block accumulator 1*1
  auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(1, 1, nbEqs));

  const std::string trsName = trs->getName();
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  std::vector<CFuint>::iterator itd;

  const std::string socketName = "boundaryMovement_" + trsName ;
  DataHandle<RealVector> boundaryDisplacement = _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();

  RealVector dirichletState(nbEqs);

  for (CFuint iState = 0; iState < trsStates->size(); ++iState) {
    State *const currState = states[(*trsStates)[iState]];

    if (currState->isParUpdatable()) {

      const CFuint localStateID = currState->getLocalID();

      if (!isUpdated[localStateID]) {

        dirichletState = boundaryDisplacement[iState];

        // coefficient i,i in global linear system
        const CFreal coeff = 1.;

        // rhs is set to
        if(_isImplicit) {
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            rhs(localStateID, iEq, nbEqs) = dirichletState[iEq] - (*currState)[iEq];
          }
        }
        else {
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            rhs(localStateID, iEq, nbEqs) = dirichletState[iEq];
          }
        }

        acc->setRowIndex(0, localStateID);

        // here we have to know how many and which vertices
        // reference the current boundary node to avoid VERY
        // expensive reallocations!!!!
        const CFuint nbEntries = bStatesNeighbors[localStateID].size();
        cf_assert(nbEntries > 0);

        for (CFuint i = 0; i < nbEntries; ++i) {
          const CFuint entryID = bStatesNeighbors[localStateID][i]->getLocalID();
          acc->setColIndex(0, entryID);
          acc->setValue(0.0);
          if (entryID == localStateID) {
            for (CFuint ib = 0; ib < nbEqs; ++ib) {
              acc->setValue(0,0,ib,ib, coeff);
            }
          }
          jacobMatrix->setValues(*acc);
        }
        isUpdated[localStateID] = true; // flagging is important!!!!!
      }
    }
  }

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();
}

//////////////////////////////////////////////////////////////////////////////

void MeshMovementDirichletBC::configure ( Config::ConfigArgs& args )
{
  FiniteElementMethodCom::configure(args);

  const std::vector<std::string>& trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];

    const std::string socketName = "boundaryMovement_" + trsName ;
    _sockets.createSocketSink<RealVector>(socketName);
  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
MeshMovementDirichletBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_states);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
