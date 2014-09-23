#include "Framework/MeshData.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/TwoLayerComputeRhs.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerComputeRhs, FluctuationSplitData, FluctSplitSpaceTimeModule> TwoLayerComputeRhsProvider("TwoLayerComputeRhs");

//////////////////////////////////////////////////////////////////////////////

TwoLayerComputeRhs::TwoLayerComputeRhs(const std::string& name) :
  ComputeRHS(name),
  _temp(),
  socket_interRhs("interRhs"),
  socket_interUpdateCoeff("interUpdateCoeff")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerComputeRhs::~TwoLayerComputeRhs()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerComputeRhs::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = ComputeRHS::needsSockets();

  result.push_back(&socket_interRhs);
  result.push_back(&socket_interUpdateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerComputeRhs::executeOnTrs()
{
  cleanRHS();

  // reset the update flag to false for all the states
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated = false;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  //  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell(); // unused
  
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  // loop over element/cell types
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    // number of cells of this type
    const CFuint nbCellsPerType = (*elementType)[iType].getNbElems();

    // loop over cells in element type
    for (CFuint iCell = 0; iCell < nbCellsPerType; ++iCell) {

      CFLogDebugMed( "Computing iCell = " << iCell << "/" << nbCellsPerType << "\n");

      // build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity& cell = *geoBuilder->buildGE();
      vector<State*> *const states = cell.getStates();

      // get the all state vectors in this cell
      getMethodData().getDistributionData().cell   = &cell;
      getMethodData().getDistributionData().cellID = cell.getID();
      getMethodData().getDistributionData().states = states;

      _fsStrategy->computeFluctuation(_residual);

      // get the all state vectors in this cell
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

      // update the rhs (it will be copied in a PetscVector later on)
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

      //release the GeometricEntity
      geoBuilder->releaseGE();

    } // end loop over cells in Type

  } // end loop over CellTypes

  // compute dynamically the CFL
  getMethodData().getCFL()->update();

}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerComputeRhs::setup()
{
  CFLogDebugMin( "TwoLayerComputeRHS::setup() BEGIN\n");
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

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  CFLogDebugMin( "TwoLayerComputeRHS::setup() END\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD
