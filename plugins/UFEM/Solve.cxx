#include "Common/BadValueException.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/Solve.hh"
#include "UFEM/ElemAssembler.hh"

//////////////////////////////////////////////////////////////////////////////
 
using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Solve, UFEMSolverData, UFEMPlugin> aSolveProvider("StdSolve");

//////////////////////////////////////////////////////////////////////////////

Solve::Solve(const std::string& name) :
  UFEMSolverCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_interStates("interStates"),
  m_nbEqs(0),
  m_dim(0)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

Solve::~Solve()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void Solve::setup()
{
  CFAUTOTRACE;

  // first call parent method
  UFEMSolverCom::setup();

  // add here specific setup
  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim   = PhysicalModelStack::getActive()->getDim();

  return;
}

//////////////////////////////////////////////////////////////////////////////

void Solve::executeOnTrs()
{
  CFAUTOTRACE;

  // watching time
  Common::Stopwatch<Common::WallTime> stopTimer;
  stopTimer.start();

  // accessing element type data
  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();

  // local vars
  const CFreal dt  = getMethodData().getDt(); // time step
  CFLog(INFO, "Time step size= " << dt << "\n");
  const CFreal invdt  = 1./dt; // time step

  // getting storages
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  // and system matrix
  SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLinearSystemSolver()[0]->getMatrix();

  // provide access to the GeometricEntities
  SafePtr< GeometricEntityPool< StdTrsGeoBuilder > > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

//   CFint pn = getMethodData().getPrintNode();
//   DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
//   State& state = *states[pn];
//   CFLog(1,"UFEM_NODAL_INFO Node: " << pn << " States: " << state << "\n");

  const CFuint nbElemTypes = elementType->size();
  const std::string trsname = trs->getName();

  // get element assembler
  ElemAssemblerMap& elmassmap = getMethodData().getelem_assemblers();

  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {

    // values of the cell states and their IDs
    const CFuint nbStates = (*elementType)[iType].getNbStates();
    RealVector nstate(nbStates*m_nbEqs);
    CFuint * nID    = new CFuint [nbStates];
    bool   * nIPU   = new bool   [nbStates];
    RealMatrix eK(m_nbEqs*nbStates,m_nbEqs*nbStates);
    RealVector eR(m_nbEqs*nbStates);

    // blockacvcumulator
    BlockAccumulator *acc = getMethodData().getLinearSystemSolver()[0]->createBlockAccumulator(nbStates,nbStates,m_nbEqs);

    // number of cells of this type
    const CFuint nbCellsPerType = (*elementType)[iType].getNbElems();

    // get element assembler
    ElemID id ( trs->getName() , iType );
    ElemAssembler *elmass= elmassmap[id];

    // call pre
    for ( ElemAssembler::VecUFEMTerm::iterator iTerm=elmass->terms.begin(); iTerm != elmass->terms.end(); ++iTerm) {
      (*iTerm)->pre();
    }

    // loop over cells in element type
    for (CFuint iCell = 0; iCell < nbCellsPerType; ++iCell)
    {
      //std::cout << "Computing iCell = " << iCell << "/" << nbCellsPerType << "\n" << std::endl;

      // build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity& cell = *geoBuilder->buildGE();

      // clean assdata
      elmass->assdata.clean();

      // get states
      vector< State* >& states = *cell.getStates();
      for (CFuint iState=0; iState<nbStates; ++iState) {
        State& state = *states[iState];
        nID[iState] = state.getLocalID();
        nIPU[iState] = state.isParUpdatable();
        for (CFuint iEq=0; iEq<m_nbEqs; ++iEq) {
          nstate[iState*m_nbEqs + iEq] = state[iEq];
        }
      }

      // prepare elemprops
      elmass->eprops->prepare(cell);
      // compute terms
      for ( ElemAssembler::VecUFEMTerm::iterator iTerm=elmass->terms.begin(); iTerm != elmass->terms.end(); ++iTerm) {
        (*iTerm)->compute(cell,elmass->assdata);
      }

      // Crank-Nicholson time stepping
      eK = invdt* elmass->assdata.T + /* 0.5 * */ elmass->assdata.A;
      eR = elmass->assdata.A*nstate + elmass->assdata.b;

      // for all the states in this cell
      for (CFuint iState=0; iState<nbStates; ++iState)
      {
        // add the contribution to the RHS vector
        for (CFuint iEq=0; iEq<m_nbEqs; ++iEq) {
          rhs(nID[iState], iEq, m_nbEqs) -= eR[iState*m_nbEqs + iEq];
        }

        // set the index of the block corresponding to the current
        // state in the jacobian matrix
        acc->setRowColIndex(iState, nID[iState]);
      }

      // add the contribution to the jacobian matrix
      acc->setValuesM(eK);
      jacobMatrix->addValues( *acc );

      //release the GeometricEntity
      geoBuilder->releaseGE();

    } // for each cell

    // call pre
    for ( ElemAssembler::VecUFEMTerm::iterator iTerm=elmass->terms.begin(); iTerm != elmass->terms.end(); ++iTerm) {
      (*iTerm)->post();
    }

    // release memory
    deletePtrArray ( nID );
    deletePtrArray ( nIPU );

  } // for each element type

  CFreal dtmult  = getMethodData().getDtmult(); // time step multiplicator
  CFreal dtlimit  = getMethodData().getDtlimit(); // time step limit
  getMethodData().setDt(min(dt*dtmult,dtlimit)); // time step

  stopTimer.stop();

  // print relative residual
  CFLog(INFO,"Assembly time: " << stopTimer << "s\n");

}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
Solve::needsSockets()
{
  CFAUTOTRACE;
  std::vector< SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM

} // namespace COOLFluiD

