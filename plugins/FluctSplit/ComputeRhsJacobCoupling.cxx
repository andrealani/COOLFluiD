#include "FluctSplit/FluctSplit.hh"
#include "ComputeRhsJacobCoupling.hh"
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

MethodCommandProvider<ComputeRhsJacobCoupling, FluctuationSplitData, FluctSplitModule>
computeRhsJacobCouplingProvider("RhsJacobCoupling");

//////////////////////////////////////////////////////////////////////////////

ComputeRhsJacobCoupling::ComputeRhsJacobCoupling(const std::string& name) :
  ComputeRHS(name),
  _lss(),
  _acc(),
  _equations(),
  _jacobStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeRhsJacobCoupling::~ComputeRhsJacobCoupling()
{
  const CFuint accSize = _acc.size();
  for (CFuint i = 0; i < accSize; ++i) {
    const CFuint nbEqs = _acc[i].size();
    for (CFuint j = 0; j < nbEqs; ++j) {
      deletePtr(_acc[i][j]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRhsJacobCoupling::executeOnTrs()
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
    for (CFuint i = 0; i < _lss.size(); ++i) {   
      SafePtr<LSSMatrix> jacobMatrix = _lss[i]->getMatrix();
      // reset to zero all non zero entries in the jacobian
      jacobMatrix->resetToZeroEntries();
    }
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
	for (CFuint iLSS = 0; iLSS < _lss.size(); ++iLSS) {   
	  // get the apropriate block matrix accumulator for this type of cell
	  BlockAccumulator *const accType = _acc[iType][iLSS];
	  
	  // set the equation subsystem descriptor which needs:
	  // 1- number of equations in the current subsystem
	  // 2- ID corresponding to the first entry in the current subsystem
	  SafePtr<vector<CFuint> > currEqs = _equations[iLSS];
	  PhysicalModelStack::getActive()->setEquationSubSysDescriptor((*currEqs)[0],  currEqs->size(), iLSS);
	  CFLog(VERBOSE, "(" << currEqs->size() << "," <<  (*currEqs)[0] << "," << iLSS <<   ")\n");
	  
	  _jacobStrategy->computeJacobianTerm(&cell,
					      *tBackResidual,
					      accType,
					      *currEqs);
	  
	  // add the values in the jacobian matrix
	  _lss[iLSS]->getMatrix()->addValues(*accType);
	  
	  // reset to zero the entries in the block accumulator
	  accType->reset();
	}
      }
      
      fsmdata.getDistributionData().isPerturb = false;
      //  std::cout<<"jacobian computed"<<std::endl;
      //release the GeometricEntity
      geoBuilder->releaseGE();
      
      // reset the equation subsystem descriptor
      PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();
    } // end loop over cells in Type
    
    nbcell+=nbCellsPerType;
  } // end loop over CellTypes

  // for safety reset the transport properties freezing to false
  diffVar->setFreezeCoeff(false);
  
  // compute dynamically the CFL if we are not doing an unsteady computation
  // Because if it is unsteady we want that it is done in the takeStepImpl
  // This is not done in unsteady because then, we can not change the CFL in the 
  // pseudo-time iterations
 /* bool isUnsteady(false);
  if(SubSystemStatusStack::getActive()->getDT() > 0.) isUnsteady = true;
  if (!isUnsteady){
    fsmdata.getCFL()->update();
  }*/
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRhsJacobCoupling::setup()
{
  CFAUTOTRACE;

  // first call parent method
  ComputeRHS::setup();
  
  const CFuint nbLSS = getMethodData().getLinearSystemSolver().size();
  
  // set the total number of equation subsystems
  PhysicalModelStack::getActive()->setTotalNbEqSS(nbLSS);

  // acquaintance of the linear systems
  _lss.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    _lss[i] = getMethodData().getLinearSystemSolver()[i];
  }
 
  _equations.resize(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    // acquaintance of the equation IDs to solve in each LSS
    _equations[i] = _lss[i]->getEquationIDs();
    cf_assert(_equations[i]->size() == _lss[i]->getNbSysEqs());
  }
  
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  _acc.resize(nbElemTypes);
  // loop over types since it can happen to deal with an hybrid mesh
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    _acc[iType].resize(nbLSS);
    // number of states in a cell of this type
    const CFuint nbStatesInType = (*elementType)[iType].getNbStates();
    for (CFuint i = 0; i < nbLSS; ++i) {
      const CFuint nbSysEq = _lss[i]->getNbSysEqs();
      // create a block matrix accumulator for this type of cell
      _acc[iType][i] = _lss[i]->createBlockAccumulator(nbStatesInType,
						       nbStatesInType,
						       nbSysEq);
    }
  }
  
  vector<vector<CFuint> > eqVarPatterns(nbLSS);
  for (CFuint i = 0; i < nbLSS; ++i) {
    eqVarPatterns[i] = *_equations[i];
  }
  
  PhysicalModelStack::getActive()->setEquationVarPatterns(eqVarPatterns);
  
  cf_assert(_fsStrategy.isNotNull());

  _jacobStrategy = getMethodData().getJacobStrategy();
  cf_assert(_jacobStrategy.isNotNull());
  
  _jacobStrategy->setup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD
