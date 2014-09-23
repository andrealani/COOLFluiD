#include "FluctSplit/FluctSplit.hh"
#include "ComputeRhsJacob.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Framework/JacobianLinearizer.hh"
#include "Framework/CFL.hh"
#include "ComputeJacobStrategy.hh"
#include "ComputeDiffusiveTerm.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/ConvergenceStatus.hh"

#include "Framework/ConvergenceMethod.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeRhsJacob, FluctuationSplitData, FluctSplitModule>
computeRhsJacobProvider("RhsJacob");

//////////////////////////////////////////////////////////////////////////////

ComputeRhsJacob::ComputeRhsJacob(const std::string& name) :
  ComputeRHS(name),
  _lss(CFNULL),
  _acc(),
  _jacobStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeRhsJacob::~ComputeRhsJacob()
{
  const CFuint accSize = _acc.size();
  for (CFuint i = 0; i < accSize; ++i) {
    deletePtr(_acc[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRhsJacob::executeOnTrs()
{
  CFAUTOTRACE;

  cleanRHS();

  // reset the update flag to false for all the states
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated = false;
  
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  if (getMethodData().doComputeJacobian()) {
    SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();
    // reset to zero all non zero entries in the jacobian
    jacobMatrix->resetToZeroEntries();
  }
  
  const CFuint nbElemTypes = elementType->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool updateDiffCoeff = true;

    FluctuationSplitData& fsmdata = getMethodData();
  DistributionData& distdata = fsmdata.getDistributionData();
  SafePtr<ConvergenceMethod> cvmth = getMethodData().getCollaborator<ConvergenceMethod>();
  SafePtr<Framework::ConvergenceMethodData>  cvmthdata = cvmth->getConvergenceMethodData();
  distdata.subiter =  cvmthdata->getConvergenceStatus().subiter;
  SafePtr<DiffusiveVarSet> diffVar = fsmdata.getDiffusiveVar();
  CFuint nbcell=0;
  // loop over element/cell types
  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    // number of cells of this type
    const CFuint nbCellsPerType = (*elementType)[iType].getNbElems();
    // get the apropriate block matrix accumulator for this type of cell
    BlockAccumulator *const accType = _acc[iType];

    Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = fsmdata.getStdTrsGeoBuilder();

    StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
    geoData.trs = getCurrentTRS();
     
    
    // loop over cells in element type
    for (CFuint iCell = 0; iCell < nbCellsPerType; ++iCell)
    {
      CFLogDebugMed( "Computing iCell = " << iCell << "/" << nbCellsPerType << "\n");

      // build the GeometricEntity
      geoData.idx = iCell+nbcell;
      GeometricEntity& cell = *geoBuilder->buildGE();
      vector<State*> *const states = cell.getStates();

      cf_assert(cell.getID() == iCell);

      distdata.cell   = &cell;
      distdata.cellID = cell.getID();
      distdata.states = states;

      const CFuint nbStatesInCell = states->size();
      // get the all state vectors in this cell
      _fsStrategy->computeFluctuation(_residual);
      //   std::cout<<"computed the advective residua;"<<std::endl;

      // transform back the residual to solution variables
      vector<RealVector> *const tBackResidual =
        _distToSolutionMatTrans->transformFromRef(&_residual);

      // for (CFuint in =0 ; in < tBackResidual->size(); ++in) {
      // 	cout << in << " node => " << (*tBackResidual)[in] << endl; 
      // }
      // abort();
      
      if (_hasArtDiffusiveTerm)
      {
	_adStrategy->addArtificialDiff(_artdiffResidual);

        for (CFuint i = 0; i < nbStatesInCell; ++i)
        {
          (*tBackResidual)[i] += _artdiffResidual[i];
        }
      }

      if (_hasDiffusiveTerm)
      {
	diffVar->setFreezeCoeff(false);

        _diffTermComputer->computeDiffusiveTerm(&cell, _diffResidual, updateDiffCoeff);
        for (CFuint i = 0; i < nbStatesInCell; ++i)
        {
          (*tBackResidual)[i] -= _diffResidual[i];
        }
      }
      //       std::cout<<"computed the full residua;"<<std::endl;
      // update the rhs (it will be copied in a Vector later on)
      // WATCH OUT: rhs has sign opposite of the one of the calculated residual
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
      {
        const CFuint stateID = (*states)[iState]->getLocalID();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
        {
          rhs(stateID, iEq, nbEqs) -= (*tBackResidual)[iState][iEq];
        }
      }
      //std::cout<<"rhs assembled"<<std::endl;

      // compute the jacobian term
      fsmdata.getDistributionData().isPerturb = true;

      diffVar->setFreezeCoeff(_freezeDiffCoeff);
      
      if (getMethodData().doComputeJacobian()) {
        _jacobStrategy->computeJacobianTerm(&cell,
                                            *tBackResidual,
                                            accType,
					    *_lss->getEquationIDs());
	
        // add the values in the jacobian matrix
        _lss->getMatrix()->addValues(*accType);
	
        // reset to zero the entries in the block accumulator
        accType->reset();
      }
      
      fsmdata.getDistributionData().isPerturb = false;
      //  std::cout<<"jacobian computed"<<std::endl;
      //release the GeometricEntity
      geoBuilder->releaseGE();

    } // end loop over cells in Type
    nbcell+=nbCellsPerType;
  } // end loop over CellTypes

  // for safety reset the transport properties freezing to false
  diffVar->setFreezeCoeff(false);
  
  // compute dynamically the CFL if we are not doing an unsteady computation
  // Because if it is unsteady we want that it is done in the takeStepImpl
  // This is not done in unsteady because then, we can not change the CFL in the 
  // pseudo-time iterations
  bool isUnsteady(false);
  if(SubSystemStatusStack::getActive()->getDT() > 0.) isUnsteady = true;
  if (!isUnsteady){
    fsmdata.getCFL()->update();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRhsJacob::setup()
{
  CFAUTOTRACE;

  // first call parent method
  ComputeRHS::setup();

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
                                               nbEqs);
  }
  
  cf_assert(_fsStrategy.isNotNull());

  _jacobStrategy = getMethodData().getJacobStrategy();
  cf_assert(_jacobStrategy.isNotNull());

  _jacobStrategy->setup();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD
