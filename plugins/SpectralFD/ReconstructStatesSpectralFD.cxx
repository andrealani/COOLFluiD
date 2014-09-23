#include "Framework/MethodStrategyProvider.hh"
#include "Framework/State.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/ReconstructStatesSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    ReconstructStatesSpectralFD,SpectralFDMethodData,ReconstructStatesSpectralFD,SpectralFDModule >
  ReconstructStatesSpectralFDProvider("ReconstructStatesSpectralFD");

//////////////////////////////////////////////////////////////////////////////

ReconstructStatesSpectralFD::ReconstructStatesSpectralFD(const std::string& name) :
  SpectralFDMethodStrategy(name),
  m_allGradMappedCoord(),
  m_gradMappedCoord()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

ReconstructStatesSpectralFD::~ReconstructStatesSpectralFD()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::reconstructState(const vector< State* >& cellStates,
                                                   State& recState,
                                                   const RealVector& recCoefs,
                                                   const vector< CFuint >& cellStateIdxsForRec)
{
  // number of cell states involved in reconstruction
  const CFuint nbrCellStatesInRec = cellStateIdxsForRec.size();
  cf_assert(nbrCellStatesInRec == recCoefs.size());
  cf_assert(nbrCellStatesInRec <= cellStates.size());

  // reconstruct the state
  recState = 0.0;
  for (CFuint iCellState = 0; iCellState < nbrCellStatesInRec; ++iCellState)
  {
    const CFuint cellStateIdx = cellStateIdxsForRec[iCellState];
    recState += recCoefs[iCellState]*(*cellStates[cellStateIdx]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::reconstructStates(const vector< State* >& cellStates,
                                                    vector< State* >& recStates,
                                                    const RealMatrix& recCoefs,
                                                    const vector< CFuint >& recStateIdxs,
                                                    const vector< CFuint >& recStateMatrixIdxs,
                                                    const vector< vector< CFuint > >& cellStateIdxs)
{
  // number of states to be reconstructed
  const CFuint nbrRecStates = recStateIdxs.size();
  cf_assert(nbrRecStates <= recStates.size());

// cout << cellStates.size() << " | " << cellStates[0]->getLocalID() << "\n" << flush;

  // reconstruct states
  for (CFuint iRecState = 0; iRecState < nbrRecStates; ++iRecState)
  {
    // index of the state that has to be reconstructed
    const CFuint recStateIdx = recStateIdxs[iRecState];

    // index in the matrix corresponding to this reconstructed state
    const CFuint recStateMatrixIdx = recStateMatrixIdxs[recStateIdx];

    // indexes of cell states involved in reconstruction
    const vector< CFuint >& cellStateIdxsForRec = cellStateIdxs[recStateIdx];

    // number of cell states involved in reconstruction
    const CFuint nbrCellStatesInRec = cellStateIdxsForRec.size();
    cf_assert(nbrCellStatesInRec <= cellStates.size());

    // dereference state to be reconstructed
    State& recState = *recStates[iRecState];

    // reconstruct the state
    recState = 0.0;
    for (CFuint iCellState = 0; iCellState < nbrCellStatesInRec; ++iCellState)
    {
      const CFuint cellStateIdx = cellStateIdxsForRec[iCellState];
      recState += recCoefs(recStateMatrixIdx,iCellState)*(*cellStates[cellStateIdx]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::reconstructExtraVars(const vector< RealVector* >& cellExtraVars,
                                                       vector< RealVector* >& recExtraVars,
                                                       const RealMatrix& recCoefs,
                                                       const vector< CFuint >& recVarIdxs,
                                                       const vector< CFuint >& recVarMatrixIdxs,
                                                       const vector< std::vector< CFuint > >& cellVarIdxs)
{
    // number of extra variables to be reconstructed
  const CFuint nbrRecVars = recVarIdxs.size();
  cf_assert(nbrRecVars <= recExtraVars.size());

  // reconstruct extra variables
  for (CFuint iRecVar = 0; iRecVar < nbrRecVars; ++iRecVar)
  {
    // index of the variable that has to be reconstructed
    const CFuint recVarIdx = recVarIdxs[iRecVar];

    // index in the matrix corresponding to this reconstructed variable
    const CFuint recVarMatrixIdx = recVarMatrixIdxs[recVarIdx];

    // indexes of cell variables involved in reconstruction
    const vector< CFuint >& cellVarIdxsForRec = cellVarIdxs[recVarIdx];

    // number of cell variables involved in reconstruction
    const CFuint nbrCellVarsInRec = cellVarIdxsForRec.size();
    cf_assert(nbrCellVarsInRec <= cellExtraVars.size());

    // dereference variable to be reconstructed
    RealVector& recVar = *recExtraVars[iRecVar];

    // reconstruct the variable
    recVar = 0.0;
    for (CFuint iCellVar = 0; iCellVar < nbrCellVarsInRec; ++iCellVar)
    {
      const CFuint cellVarIdx = cellVarIdxsForRec[iCellVar];
      cf_assert(recVar.size() == cellExtraVars[cellVarIdx]->size());
      recVar += recCoefs(recVarMatrixIdx,iCellVar)*(*cellExtraVars[cellVarIdx]);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::reconstructGradients(const vector< vector< RealVector >* >& cellGradients,
                                                       vector< vector< RealVector* > >& recGradients,
                                                       const RealMatrix& recCoefs,
                                                       const vector< CFuint >& recGradIdxs,
                                                       const vector< CFuint >& recGradMatrixIdxs,
                                                       const vector< vector< CFuint > >& cellGradIdxs)
{
  // number of gradients to be reconstructed
  const CFuint nbrRecGrads = recGradIdxs.size();
  cf_assert(nbrRecGrads <= recGradients.size());

  // number of gradient variables
  cf_assert(cellGradients.size() > 0);
  const CFuint nbrGradVars = cellGradients[0]->size();

  // reconstruct gradients
  for (CFuint iRecGrad = 0; iRecGrad < nbrRecGrads; ++iRecGrad)
  {
    // index of the gradient that has to be reconstructed
    const CFuint recGradIdx = recGradIdxs[iRecGrad];

    // index in the matrix corresponding to this reconstructed gradient
    const CFuint recGradMatrixIdx = recGradMatrixIdxs[recGradIdx];

    // indexes of cell gradients involved in reconstruction
    const vector< CFuint >& cellGradIdxsForRec = cellGradIdxs[recGradIdx];

    // number of cell gradients involved in reconstruction
    const CFuint nbrCellGradsInRec = cellGradIdxsForRec.size();
    cf_assert(nbrCellGradsInRec <= cellGradients.size());

    // loop over the different variables in the gradient
    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      // dereference gradient to be reconstructed
      RealVector& recGrad = *recGradients[iRecGrad][iGradVar];

      // reconstruct the gradient
      recGrad = 0.0;
      for (CFuint iCellGrad = 0; iCellGrad < nbrCellGradsInRec; ++iCellGrad)
      {
        const CFuint cellGradIdx = cellGradIdxsForRec[iCellGrad];
        recGrad += recCoefs(recGradMatrixIdx,iCellGrad)*(*cellGradients[cellGradIdx])[iGradVar];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::reconstructPhysVar(const CFuint iVar,
                                                     const vector< State* >& cellStates,
                                                     vector< State* >& recStates,
                                                     const RealMatrix& recCoefs,
                                                     const vector< CFuint >& recStateIdxs,
                                                     const vector< CFuint >& recStateMatrixIdxs,
                                                     const vector< vector< CFuint > >& cellStateIdxs)
{
  // number of states to be reconstructed
  const CFuint nbrRecStates = recStateIdxs.size();
  cf_assert(nbrRecStates <= recStates.size());

  // reconstruct states
  for (CFuint iRecState = 0; iRecState < nbrRecStates; ++iRecState)
  {
    // index of the state that has to be reconstructed
    const CFuint recStateIdx = recStateIdxs[iRecState];

    // index in the matrix corresponding to this reconstructed state
    const CFuint recStateMatrixIdx = recStateMatrixIdxs[recStateIdx];

    // indexes of cell states involved in reconstruction
    const vector< CFuint >& cellStateIdxsForRec = cellStateIdxs[recStateIdx];

    // number of cell states involved in reconstruction
    const CFuint nbrCellStatesInRec = cellStateIdxsForRec.size();
    cf_assert(nbrCellStatesInRec <= cellStates.size());

    // dereference state to be reconstructed
    State& recState = *recStates[iRecState];

    // reconstruct the variable
    recState[iVar] = 0.0;
    for (CFuint iCellState = 0; iCellState < nbrCellStatesInRec; ++iCellState)
    {
      const CFuint cellStateIdx = cellStateIdxsForRec[iCellState];
      recState[iVar] += recCoefs(recStateMatrixIdx,iCellState)*(*cellStates[cellStateIdx])[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::reconstructPhysVarGrad(const CFuint iVar, const vector< vector< RealVector >* >& cellGradients,
                                                         vector< vector< RealVector* > >& recGradients,
                                                         const RealMatrix& recCoefs,
                                                         const vector< CFuint >& recGradIdxs,
                                                         const vector< CFuint >& recGradMatrixIdxs,
                                                         const vector< vector< CFuint > >& cellGradIdxs)
{
  // number of gradients to be reconstructed
  const CFuint nbrRecGrads = recGradIdxs.size();
  cf_assert(nbrRecGrads <= recGradients.size());
  cf_assert(cellGradients.size() > iVar);

  // reconstruct gradients
  for (CFuint iRecGrad = 0; iRecGrad < nbrRecGrads; ++iRecGrad)
  {
    // index of the gradient that has to be reconstructed
    const CFuint recGradIdx = recGradIdxs[iRecGrad];

    // index in the matrix corresponding to this reconstructed gradient
    const CFuint recGradMatrixIdx = recGradMatrixIdxs[recGradIdx];

    // indexes of cell gradients involved in reconstruction
    const vector< CFuint >& cellGradIdxsForRec = cellGradIdxs[recGradIdx];

    // number of cell gradients involved in reconstruction
    const CFuint nbrCellGradsInRec = cellGradIdxsForRec.size();
    cf_assert(nbrCellGradsInRec <= cellGradients.size());

    // dereference gradient to be reconstructed
    RealVector& recGrad = *recGradients[iRecGrad][iVar];

    // reconstruct the gradient
    recGrad = 0.0;
    for (CFuint iCellGrad = 0; iCellGrad < nbrCellGradsInRec; ++iCellGrad)
    {
      const CFuint cellGradIdx = cellGradIdxsForRec[iCellGrad];
      recGrad += recCoefs(recGradMatrixIdx,iCellGrad)*(*cellGradients[cellGradIdx])[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::computePolyGradients(const vector< Framework::State* >& cellStates,
                                                       vector< vector< RealVector* > >& recGradients,
                                                       const vector< vector< vector< CFreal > > >& derivCoefs,
                                                       const vector< CFuint >& recGradIdxs,
                                                       const vector< vector< vector< CFuint > > >& cellStateIdxs,
                                                       const vector< RealMatrix >& invJacobMatr)
{
  // number of gradients to be reconstructed
  const CFuint nbrRecGrads = recGradIdxs.size();
  cf_assert(nbrRecGrads <= recGradients.size());

  // number of physical variables and dimensionality
  cf_assert(cellStates.size() > 0);
  const CFuint nbrPhysVars = cellStates[0]->size();
  cf_assert(recGradients.size() > 0);
  cf_assert(recGradients[0].size() > 0);
  const CFuint dim = recGradients[0][0]->size();

  // reconstruct gradients
  for (CFuint iRecGrad = 0; iRecGrad < nbrRecGrads; ++iRecGrad)
  {
    // dereference gradient to be computed
    vector< RealVector* >& currRecGrad = recGradients[iRecGrad];

    // index of the gradient that has to be reconstructed
    const CFuint recGradIdx = recGradIdxs[iRecGrad];

    // coefficients for gradient computation
    const vector< vector< CFreal > >& currDerivCoefs = derivCoefs[recGradIdx];
    cf_assert(currDerivCoefs.size() == dim);

    // indexes of corresponding cell states
    const vector< vector< CFuint > >& currDerivIdxs = cellStateIdxs[recGradIdx];
    cf_assert(currDerivIdxs.size() == dim);

    // reconstruct the gradient with respect to the mapped coordinates
    for (CFuint iDir = 0; iDir < dim; ++iDir)
    {
      // coefficients and solution point indexes for current direction
      const vector< CFreal >& dirDerivCoefs = currDerivCoefs[iDir];
      const vector< CFuint >& dirDerivIdxs  = currDerivIdxs [iDir];
      cf_assert(dirDerivCoefs.size() == dirDerivIdxs.size());

      // number of cell states involved
      const CFuint nbrCellStates = dirDerivCoefs.size();

      // loop over the different variables in the gradient
      for (CFuint iVar = 0; iVar < nbrPhysVars; ++iVar)
      {
        // dereference gradient to be reconstructed
        RealVector& recGrad = m_allGradMappedCoord[iVar];
        recGrad[iDir] = 0.0;
        for (CFuint iCellState = 0; iCellState < nbrCellStates; ++iCellState)
        {
          const CFuint stateIdx = dirDerivIdxs[iCellState];
          recGrad[iDir] += dirDerivCoefs[iCellState]*(*cellStates[stateIdx])[iVar];
        }
      }
    }

    // transform to global coordinates
    const RealMatrix& locInvJacobMatr = invJacobMatr[iRecGrad];
    for (CFuint iVar = 0; iVar < nbrPhysVars; ++iVar)
    {
      *currRecGrad[iVar] = locInvJacobMatr*m_allGradMappedCoord[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::computePolyGradient(const CFuint iVar,
                                                      const vector< Framework::State* >& cellStates,
                                                      vector< vector< RealVector* > >& recGradients,
                                                      const vector< vector< vector< CFreal > > >& derivCoefs,
                                                      const vector< CFuint >& recGradIdxs,
                                                      const vector< vector< vector< CFuint > > >& cellStateIdxs,
                                                      const vector< RealMatrix >& invJacobMatr)
{
  // number of gradients to be reconstructed
  const CFuint nbrRecGrads = recGradIdxs.size();
  cf_assert(nbrRecGrads <= recGradients.size());

  // number of physical variables and dimensionality
  cf_assert(cellStates.size() > 0);
  cf_assert(cellStates[0]->size() > iVar);
  cf_assert(recGradients.size() > 0);
  cf_assert(recGradients[0].size() > 0);
  const CFuint dim = recGradients[0][0]->size();

  // reconstruct gradients
  for (CFuint iRecGrad = 0; iRecGrad < nbrRecGrads; ++iRecGrad)
  {
    // dereference gradient to be computed
    vector< RealVector* >& currRecGrad = recGradients[iRecGrad];

    // index of the gradient that has to be reconstructed
    const CFuint recGradIdx = recGradIdxs[iRecGrad];

    // coefficients for gradient computation
    const vector< vector< CFreal > >& currDerivCoefs = derivCoefs[recGradIdx];
    cf_assert(currDerivCoefs.size() == dim);

    // indexes of corresponding cell states
    const vector< vector< CFuint > >& currDerivIdxs = cellStateIdxs[recGradIdx];
    cf_assert(currDerivIdxs.size() == dim);

    // reconstruct the gradient
    for (CFuint iDir = 0; iDir < dim; ++iDir)
    {
      // coefficients and solution point indexes for current direction
      const vector< CFreal >& dirDerivCoefs = currDerivCoefs[iDir];
      const vector< CFuint >& dirDerivIdxs  = currDerivIdxs [iDir];
      cf_assert(dirDerivCoefs.size() == dirDerivIdxs.size());

      // number of cell states involved
      const CFuint nbrCellStates = dirDerivCoefs.size();

      // compute gradient
      m_gradMappedCoord[iDir] = 0.0;
      for (CFuint iCellState = 0; iCellState < nbrCellStates; ++iCellState)
      {
        const CFuint stateIdx = dirDerivIdxs[iCellState];
        m_gradMappedCoord[iDir] += dirDerivCoefs[iCellState]*(*cellStates[stateIdx])[iVar];
      }
    }

    // transform to global coordinates
    *currRecGrad[iVar] = invJacobMatr[iRecGrad]*m_gradMappedCoord;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFD::setup()
{
  SpectralFDMethodStrategy::setup();

  const CFreal dim   = PhysicalModelStack::getActive()->getDim ();
  const CFreal nbEqs = PhysicalModelStack::getActive()->getNbEq();

  m_gradMappedCoord.resize(dim);
  m_allGradMappedCoord.resize(nbEqs,RealVector(dim));
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
