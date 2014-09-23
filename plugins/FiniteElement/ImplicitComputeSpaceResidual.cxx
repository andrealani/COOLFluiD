#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"

#include "FiniteElement/ComputeResidualStrategy.hh"
#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/ImplicitComputeSpaceResidual.hh"

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

MethodCommandProvider<ImplicitComputeSpaceResidual, FiniteElementMethodData, FiniteElementModule> implicitComputeSpaceResidualProvider("ImplicitComputeSpaceResCom");

//////////////////////////////////////////////////////////////////////////////

ImplicitComputeSpaceResidual::ImplicitComputeSpaceResidual
(const std::string& name) :
  ComputeSpaceResidual(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ImplicitComputeSpaceResidual::~ImplicitComputeSpaceResidual()
{
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidual::setup()
{
  CFAUTOTRACE;

  // first call parent method
  ComputeSpaceResidual::setup();

  // add here specific setup
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidual::executeOnTrs()
{
  CFAUTOTRACE;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  FiniteElementMethodData& femdata = getMethodData();

  SafePtr<LSSMatrix> jacobMatrix =
    femdata.getLinearSystemSolver()[0]->getMatrix();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = femdata.getStdTrsGeoBuilder();

  Common::SafePtr<ComputeJacobStrategy> jacob_strategy =
    femdata.getJacobianStrategy();

  Common::SafePtr<ComputeResidualStrategy> rhs_strategy =
    femdata.getResidualStrategy();

  LocalElementData& local_elem_data = femdata.getLocalElementData();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData(getCurrentTRS()->getName());

  const CFuint nbElemTypes    = elementType->size();
  const CFuint nbEqs          = PhysicalModelStack::getActive()->getNbEq();

  // loop over element/cell types
  local_elem_data.trs     = getCurrentTRS();
  PhysicalModelStack::getActive()->setCurrentZone(getCurrentTRS()->getName());

  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

    const CFuint nbStatesInCell = (*elementType)[iType].getNbStates();

    FElemTypeData elemTData = m_map_femdata.find(nbStatesInCell);

    local_elem_data.nbEqs     = nbEqs;
    local_elem_data.nbStates  = nbStatesInCell;
    local_elem_data.blockacc  = elemTData.first;
    local_elem_data.stiff_mat = elemTData.second;
    local_elem_data.load_vec  = elemTData.third;
    local_elem_data.residual  = elemTData.fourth;

    BlockAccumulator&   acc      = *local_elem_data.blockacc;
    vector<RealVector>& residual = *local_elem_data.residual;

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
      acc.reset();

      rhs_strategy->computeElementResidual(residual);
      // write to the right hand side
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
      {
        const CFuint stateID = states[iState]->getLocalID();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          rhs(stateID, iEq, nbEqs) -= residual[iState][iEq];
        }
      }

      // compute the jacobian matrix if it is not frozen
      if(!femdata.isSysMatrixFrozen())
      {
        jacob_strategy->computeJacobianTerm();
        jacobMatrix->addValues(acc); // add the values to the jacoban matrix
      }
      geoBuilder->releaseGE(); // release the GeometricEntity
    } // end loop cells
  } // end loop element types
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
