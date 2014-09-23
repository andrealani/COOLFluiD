#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/ExplicitComputeSpaceResidual.hh"
#include "FiniteElement/ComputeResidualStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ExplicitComputeSpaceResidual, FiniteElementMethodData, FiniteElementModule> explicitComputeSpaceResidualProvider("ExplicitComputeSpaceResCom");

//////////////////////////////////////////////////////////////////////////////

ExplicitComputeSpaceResidual::ExplicitComputeSpaceResidual(const std::string& name) :
  ComputeSpaceResidual(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ExplicitComputeSpaceResidual::~ExplicitComputeSpaceResidual()
{
}

//////////////////////////////////////////////////////////////////////////////

void ExplicitComputeSpaceResidual::setup()
{
  CFAUTOTRACE;

  // first call parent method
  ComputeSpaceResidual::setup();

  // add here specific setup
}

//////////////////////////////////////////////////////////////////////////////

void ExplicitComputeSpaceResidual::executeOnTrs()
{
  CFAUTOTRACE;

  FiniteElementMethodData& femdata  = getMethodData();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<LinearSystemSolver> lss = femdata.getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = femdata.getStdTrsGeoBuilder();

  LocalElementData& local_elem_data = femdata.getLocalElementData();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData(trs->getName());

  const CFuint nbElemTypes    = elementType->size();
  const CFuint nbEqs          = PhysicalModelStack::getActive()->getNbEq();

  local_elem_data.trs     = getCurrentTRS();
  PhysicalModelStack::getActive()->setCurrentZone(getCurrentTRS()->getName());
  // loop over element/cell types
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

    const CFuint nbStatesInCell = (*elementType)[iType].getNbStates();

    FElemTypeData elemTData = m_map_femdata.find(nbStatesInCell);

    local_elem_data.nbEqs     = nbEqs;
    local_elem_data.nbStates  = nbStatesInCell;
    local_elem_data.blockacc  = elemTData.first;
    local_elem_data.stiff_mat = elemTData.second;
    local_elem_data.load_vec  = elemTData.third;
    local_elem_data.residual  = elemTData.fourth;

    RealVector& elemVec = *local_elem_data.load_vec;
    RealMatrix& elemMat = *local_elem_data.stiff_mat;

    BlockAccumulator&   acc      = *local_elem_data.blockacc;
// unused //    vector<RealVector>& residual = *local_elem_data.residual;

    // get number of cells of this type
    const CFuint nbCellsPerType = (*elementType)[iType].getNbElems();
    // loop over cells in element type
    for (CFuint iCell = 0; iCell < nbCellsPerType; ++iCell) {

      CFLogDebugMed( "Computing iCell = " << iCell << "/" << nbCellsPerType << "\n");

      // build the GeometricEntity
      geoData.idx = iCell;
      local_elem_data.cell  = geoBuilder->buildGE();
      GeometricEntity& cell = *local_elem_data.cell;

      getMethodData().getFEMVolumeIntegrator()->precomputeCellData(&cell);

      vector<State*>& states = *cell.getStates();
      cf_assert( nbStatesInCell == states.size() );

      // compute rhs
      femdata.getResidualStrategy()->computeElemVector();

      for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
      {
        const CFuint stateID = states[iState]->getLocalID();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) += elemVec[iState*nbEqs + iEq];
        }
      }

      // compute the system matrix if it is not frozen
      if(!femdata.isSysMatrixFrozen())
      {
        // set the IDs on the blockaccumulator
        for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
        {
          const CFuint stateID = states[iState]->getLocalID();
          acc.setRowColIndex(iState, stateID);
        }

        // compute element matrix and vector
        femdata.getResidualStrategy()->computeElemMatrix();

        acc.setValuesM(elemMat);

        // add the values in the jacobian matrix
        jacobMatrix->addValues(acc);

      } // frozen?

      geoBuilder->releaseGE(); // release the GeometricEntity

    } // end loop cells

  } // end loop element types

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

#if 0
      // Compute and  Insert result into the Block Accumulator
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {

        if ((*cellStates)[iState]->isParUpdatable()) {

        State *const currState = (*cellStates)[iState];
        CFuint stateID = currState->getLocalID();
        //CFout << "iState: " << iState << "\n";
        // set the index of the block corresponding to the current
        // state in the jacobian matrix
        acc->setRowColIndex(iState, currState->getLocalID());

        for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
          //CFout << "jState: " << iState << "\n";

            // Convective Term

            // Diffusive Term
            //CFout << "Start Diffusive Term Computation" << "\n";
            diffTermComputer->setIndexes(iState,jState);
            diffTermComputer->computeTerm(cell, _integResult);

            acc->addValues(iState,
                           jState,
                           _integResult);
            //CFout << "End Diffusive Term Computation" << "\n";
            // Linear Source Term
            //CFout << "Start Linear Source Term Computation" << "\n";
            linearSourceTermComputer->setIndexes(iState,jState);
            linearSourceTermComputer->computeTerm(cell, _integResult);

            CFLogDebugMax("Adding values to " << iState << " " << jState << "\n");

            acc->addValues(iState,
                           jState,
                           _integResult);
            //CFout << "End Linear Source Term Computation" << "\n";
          }

        // Independent Source Term
        //CFout << "Start Indep Source Term Computation" << "\n";
        indepSourceTermComputer->setIndexes(iState,0);
        indepSourceTermComputer->computeTerm(cell, _integResultVec);
        // the independant source term is added to the rhs
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
           _rhs(stateID, iEq, nbEqs) += _integResultVec[iEq];
        }
        //CFout << "End Indep Source Term Computation" << "\n";
        }
      }
#endif
