#include "Common/OSystem.hh"
#include "Common/ProcessInfo.hh"

#include "Framework/MeshData.hh"
#include "Framework/JacobianLinearizer.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/ArtificialDiffusionStrategy.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/ComputeRHS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ComputeRHS,
                       FluctuationSplitData,
                       FluctSplitModule>
aStdComputeRHSProvider("StdComputeRHS");

//////////////////////////////////////////////////////////////////////////////

ComputeRHS::ComputeRHS(const std::string& name) :
  FluctuationSplitCom(name),
  _solutionToDistMatTrans(CFNULL),
  _distToSolutionMatTrans(CFNULL),
  _linearToDistMatTrans(CFNULL),
  _solutionToLinearInUpdateMatTrans(CFNULL),
  _solutionToLinearMatTrans(CFNULL),
  _solutionToUpdateMatTrans(CFNULL),
  _updateToLinearVecTrans(CFNULL),
  _distribVar(CFNULL),
  _updateVar(CFNULL),
  _linearVar(CFNULL),
  _diffTermComputer(CFNULL),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_normals("normals"),
  socket_isBState("isBState"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  _residual(0),
  _diffResidual(0),
  _artdiffResidual(0),
  _fsStrategy(CFNULL),
  _adStrategy(CFNULL),
  _hasDiffusiveTerm(false),
  _hasArtDiffusiveTerm(false)
{
  addConfigOptionsTo(this);

  _freezeDiffCoeff = false;
  setParameter("FreezeDiffCoeff",&_freezeDiffCoeff);
}

//////////////////////////////////////////////////////////////////////////////

ComputeRHS::~ComputeRHS()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRHS::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeRHS::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_normals);
  result.push_back(&socket_isBState);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRHS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool > ("FreezeDiffCoeff", "Flag forcing to freeze diffusive coefficients");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRHS::setup()
{
  CFAUTOTRACE;
  FluctuationSplitCom::setup();

  // get the distribution var set
  _distribVar = getMethodData().getDistribVar();
  // get the update var set
  _updateVar = getMethodData().getUpdateVar();
  // get the linearization var set
  _linearVar = getMethodData().getLinearVar();

  _diffTermComputer = getMethodData().getDiffusiveTermComputer();

  _solutionToDistMatTrans = getMethodData().getSolutionToDistribMatTrans();
  _distToSolutionMatTrans = getMethodData().getDistribToSolutionMatTrans();
  _linearToDistMatTrans = getMethodData().getLinearToDistribMatTrans();
  _solutionToLinearInUpdateMatTrans = getMethodData().getSolutionToLinearInUpdateMatTrans();
  _solutionToLinearMatTrans = getMethodData().getSolutionToLinearMatTrans();
  _updateToLinearVecTrans = getMethodData().getUpdateToLinearVecTrans();
  _solutionToUpdateMatTrans = getMethodData().getSolToUpdateInUpdateMatTrans();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // allocating data for the temporary local residual
  _residual.resize(maxNbStatesInCell);
  for (CFuint i = 0; i < maxNbStatesInCell; ++i)
  {
    _residual[i].resize(nbEqs);
    _residual[i] = 0.0;
  }

  // allocating data for the diffusive residual
  _diffResidual.resize(maxNbStatesInCell);
  for (CFuint i = 0; i < maxNbStatesInCell; ++i)
  {
    _diffResidual[i].resize(nbEqs);
    _diffResidual[i] = 0.0;
  }

  // allocating data for the artificial dissipation residual
  _artdiffResidual.resize(maxNbStatesInCell);
  for (CFuint i = 0; i < maxNbStatesInCell; ++i)
  {
    _artdiffResidual[i].resize(nbEqs);
    _artdiffResidual[i] = 0.0;
  }

  // access to the FS strategy
  _fsStrategy = getMethodData().getFluctSplitStrategy();
  cf_assert(_fsStrategy.isNotNull());

  // setup artificial diffusion strategy
  _adStrategy = getMethodData().getArtificialDiffusionStrategy();
  cf_assert(_adStrategy.isNotNull());

  // flag telling if a diffusive term has to be computed
  _hasDiffusiveTerm = !(_diffTermComputer->isNull());

  // flag telling if a artificial diffusive term has to be computed
  _hasArtDiffusiveTerm = !(_adStrategy->isNull());
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRHS::cleanRHS()
{
  Framework::DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  rhs = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRHS::executeOnTrs()
{
  CFAUTOTRACE;

  cf_assert(isSetup());

  // reset the RHS to 0.
  cleanRHS();
  // reset the update flag to 0 for all the states
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  isUpdated = false;

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool updateDiffCoeff = true;
  DataHandle<State*,GLOBAL> ss = socket_states.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = getCurrentTRS();

  // prepares to loop over cells by getting the GeometricEntityPool
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();
  const CFuint nbGeos = cells->getLocalNbGeoEnts();

  DistributionData& ddata = getMethodData().getDistributionData();
  cf_assert( _fsStrategy.isNotNull() ); // certify that the split strategy is not null

  // LOOP on Geometric Entities (elements)
  for (CFuint cellID = 0; cellID < nbGeos; ++cellID)
  {
//     CFout << "+++  Cell [" << cellID << "]\n";

    // build the GeometricEntity
    geoData.idx = cellID;
    GeometricEntity& cell = *geoBuilder->buildGE();
    vector<State*> *const states = cell.getStates();
    const CFuint nbStatesInCell = states->size();
    
    ddata.cell   = &cell;
    ddata.cellID = cell.getID();
    ddata.states = states;
    cf_assert(cell.getID() == cellID);
        
    _fsStrategy->computeFluctuation(_residual);

//     for ( CFuint i = 0; i < _residual.size() ; ++i )
//     { CF_DEBUG_OBJ ( _residual[i] ); }
//     CF_DEBUG_EXIT;

    // transform the residual back from the
    // distribution variables to the solution variables
    vector<RealVector> *const tBackResidual =
      _distToSolutionMatTrans->transformFromRef(&_residual);

    /// @todo fireup the diffusive terms (includes both artificial and physical terms)


    if (_hasArtDiffusiveTerm)
    {
      _adStrategy->setCell(&cell);
      _adStrategy->addArtificialDiff(_artdiffResidual);

      for (CFuint i = 0; i < nbStatesInCell; ++i)
      {
        (*tBackResidual)[i] += _artdiffResidual[i];
      }
    }

    if (_hasDiffusiveTerm)
    {
      _diffTermComputer->computeDiffusiveTerm(&cell,
                                              _diffResidual,
                                              updateDiffCoeff);

      for (CFuint i = 0; i < nbStatesInCell; ++i)
      {
        (*tBackResidual)[i] -= _diffResidual[i];
      }
    }

    for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
    {
      const CFuint stateID = (*states)[iState]->getLocalID();
      CFLogDebugMax( "rhs = ");

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        rhs(stateID, iEq, nbEqs) -= (*tBackResidual)[iState][iEq];
        CFLogDebugMax( rhs(stateID, iEq, nbEqs) << " ");
      }
      CFLogDebugMax( "\n");
    }


    //release the GeometricEntity
    geoBuilder->releaseGE();
  }

  // transform the residual from the solution variables
  // to the update variables if needed
  if (getMethodData().isResidualTransformationNeeded()) {
    transformResidual();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeRHS::transformResidual()
{
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();
  RealVector tempRes(nbEqs);
  RealVector res(nbEqs);

  for(CFuint iState = 0; iState < nbStates; ++iState) {
    // set and get the transformation matrix in the update variables
    _solutionToUpdateMatTrans->setMatrix(*states[iState]);
    const RealMatrix& matrix = *_solutionToUpdateMatTrans->getMatrix();

    // copy the rhs in a given temporary array
    const CFuint startID = iState*nbEqs;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      tempRes[iEq] = rhs[startID + iEq];
    }

    // compute the transformed residual
    res = matrix*tempRes;

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      rhs[startID + iEq] = res[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD
