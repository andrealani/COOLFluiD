#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "TwoLayerComputeRhsJacobImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Framework/JacobianLinearizer.hh"
#include "Framework/CFL.hh"
#include "ComputeJacobStrategy.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerComputeRhsJacobImpl, FluctuationSplitData, FluctSplitSpaceTimeModule> TwoLayerComputeRhsJacobImplProvider("TwoLayerRhsJacobImpl");

//////////////////////////////////////////////////////////////////////////////

TwoLayerComputeRhsJacobImpl::TwoLayerComputeRhsJacobImpl(const std::string& name) :
  ComputeRHS(name),
  _lss(CFNULL),
  _jacobStrategy(),
  _acc(),
  _temp(),
  socket_interRhs("interRhs"),
  socket_interUpdateCoeff("interUpdateCoeff")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerComputeRhsJacobImpl::~TwoLayerComputeRhsJacobImpl()
{
  const CFuint accSize = _acc.size();
  for (CFuint i = 0; i < accSize; ++i) {
    deletePtr(_acc[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerComputeRhsJacobImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = ComputeRHS::needsSockets();

  result.push_back(&socket_interRhs);
  result.push_back(&socket_interUpdateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerComputeRhsJacobImpl::executeOnTrs()
{
  cleanRHS();

  // reset the update flag to false for all the states
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated = false;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();
  // reset to zero all non zero entries in the jacobian
  jacobMatrix->resetToZeroEntries();

  const CFuint nbElemTypes = elementType->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // unused // const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  // loop over element/cell types
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    CFLogDebugMed( "Cell Type = " << iType << "\n");
    // number of cells of this type
    const CFuint nbCellsPerType = (*elementType)[iType].getNbElems();
    // get the apropriate block matrix accumulator for this type of cell
    BlockAccumulator *const accType = _acc[iType];

    // loop over cells in element type
    for (CFuint iCell = 0; iCell < nbCellsPerType; ++iCell) {

      CFLogDebugMed( "Computing iCell = " << iCell << "/" << nbCellsPerType << "\n");

      // build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity& cell = *geoBuilder->buildGE();

      // get the all state vectors in this cell
      vector<State*> *const states = cell.getStates();

      getMethodData().getDistributionData().cell   = &cell;
      getMethodData().getDistributionData().cellID = cell.getID();
      cf_assert(cell.getID() == iCell);
      getMethodData().getDistributionData().states = states;

      // get the all state vectors in this cell
      _fsStrategy->computeFluctuation(_residual);

      const CFuint nbStatesInCell = states->size();

      // Pointer to the residual (intermediate)
      for (CFuint i=0; i < nbStatesInCell; ++i){
        for (CFuint j=0; j < nbEqs; ++j){
          _temp[i][j] = _residual[i][j];
        }
      }

      // transform back the residual to solution variables
      vector<RealVector> * tBackResidual =
        _distToSolutionMatTrans->transformFromRef(&_temp);

      // update the interRhs (it will be copied in a rhs LSSVector later on)
      // WATCH OUT: rhs has sign opposite of the one of the calculated residual
      CFLogDebugMax( "interRhs = " <<"\n");
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        const CFuint stateID = (*states)[iState]->getLocalID();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          interRhs(stateID, iEq, nbEqs) -= (*tBackResidual)[iState][iEq];
          CFLogDebugMax( interRhs(stateID, iEq, nbEqs) << " ");
        }
        CFLogDebugMax( "\n");
      }

      /// Pointer to the residual (current)
      for (CFuint i=0; i < nbStatesInCell; ++i){
        for (CFuint j=0; j < nbEqs; ++j){
          _temp[i][j] = _residual[i][j+nbEqs];
        }
      }

      // transform back the residual to solution variables
      tBackResidual = _distToSolutionMatTrans->transformFromRef(&_temp);

      // update the rhs (it will be copied in a LSSVector later on)
      // WATCH OUT: rhs has sign opposite of the one of the calculated residual
      CFLogDebugMax( "rhs = " <<"\n");
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        const CFuint stateID = (*states)[iState]->getLocalID();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) -= (*tBackResidual)[iState][iEq];

          CFLogDebugMax( rhs(stateID, iEq, nbEqs) << " ");
        }
        CFLogDebugMax( "\n");
      }


      getMethodData().getDistributionData().isPerturb = true;

      // compute the jacobian term
      _jacobStrategy->computeJacobianTerm(&cell,
                                          _residual,
                                          accType,
					  *_lss->getEquationIDs());
      //states, _residual, accType);
      
      // add the values in the jacobian matrix
      jacobMatrix->addValues(*accType);

      // reset to zero the entries in the block accumulator
      accType->reset();

      getMethodData().getDistributionData().isPerturb = false;

      //release the GeometricEntity
      geoBuilder->releaseGE();
    } // end loop over cells in Type

  } // end loop over CellTypes

  // compute dynamically the CFL
  getMethodData().getCFL()->update();
 
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerComputeRhsJacobImpl::setup()
{
  // first call parent method
  ComputeRHS::setup();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  /// Modifying the size of the cell residual (for two layers)
  _residual.resize(maxNbStatesInCell);
  // allocating data for the temporary local residual
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _residual[i].resize(PhysicalModelStack::getActive()->getNbEq()*2);
    _residual[i] = 0.0;
  }

  // Sets the size of _temp
  _temp.resize(maxNbStatesInCell);
  // allocating data for the temporary local residual
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _temp[i].resize(PhysicalModelStack::getActive()->getNbEq());
  }

  _lss = getMethodData().getLinearSystemSolver()[0];

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  _acc.resize(nbElemTypes);
  // loop over types since it can happen to deal with an hybrid mesh
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    // number of states in a cell of this type
    const CFuint nbStatesInType = (*elementType)[iType].getNbStates();
    // create a block matrix accumulator for this type of cell
    _acc[iType] = _lss->createBlockAccumulator(nbStatesInType,
                  nbStatesInType,
                  nbEqs*2);
  }

  cf_assert(_fsStrategy.isNotNull());

  _jacobStrategy = getMethodData().getJacobStrategy();
  cf_assert(_jacobStrategy.isNotNull());

  _jacobStrategy->setup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
