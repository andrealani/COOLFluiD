#include "Framework/MethodStrategyProvider.hh"
#include "Framework/State.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/ReconstructStatesSpectralFV.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    ReconstructStatesSpectralFV,SpectralFVMethodData,ReconstructStatesSpectralFV,SpectralFVModule >
  ReconstructStatesSpectralFVProvider("ReconstructStatesSpectralFV");

//////////////////////////////////////////////////////////////////////////////

ReconstructStatesSpectralFV::ReconstructStatesSpectralFV(const std::string& name) :
  SpectralFVMethodStrategy(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

ReconstructStatesSpectralFV::~ReconstructStatesSpectralFV()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFV::reconstructState(const vector< State* >& cellStates,
                                                   State& recState,
                                                   const vector< CFreal >& recCoefs,
                                                   const CFuint nbrCellStates)
{
  cf_assert(nbrCellStates <= cellStates.size());
  cf_assert(nbrCellStates <= recCoefs.size()  );

  // reconstruct the state
  recState = 0.0;
  for (CFuint iCellState = 0; iCellState < nbrCellStates; ++iCellState)
  {
    recState += recCoefs[iCellState]*(*cellStates[iCellState]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFV::reconstructStates(const vector< State* >& cellStates,
                                                    vector< State* >& recStates,
                                                    const vector< vector< CFreal > >& recCoefs,
                                                    const CFuint nbrCellStates)
{
  // number of states to be reconstructed
  const CFuint nbrRecStates = recCoefs.size();
  cf_assert(nbrRecStates <= recStates.size());
  cf_assert(nbrCellStates <= cellStates.size());

  // reconstruct states
  for (CFuint iRecState = 0; iRecState < nbrRecStates; ++iRecState)
  {
    cf_assert(nbrCellStates <= recCoefs[iRecState].size());

    // dereference state to be reconstructed
    State& recState = *recStates[iRecState];

    // reconstruct the state
    recState = 0.0;
    for (CFuint iCellState = 0; iCellState < nbrCellStates; ++iCellState)
    {
      recState += recCoefs[iRecState][iCellState]*(*cellStates[iCellState]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFV::reconstructExtraVars(const vector< RealVector* >& cellExtraVars,
                                                       vector< RealVector* >& recExtraVars,
                                                       const vector< std::vector< CFreal > >& recCoefs,
                                                       const CFuint nbrCellExtraVars)
{
  // number of extra variables to be reconstructed
  const CFuint nbrRecExtraVars = recCoefs.size();
  cf_assert(nbrRecExtraVars<= recExtraVars.size());
  cf_assert(nbrCellExtraVars <= cellExtraVars.size());

  // reconstruct extra variables
  for (CFuint iRecVar = 0; iRecVar < nbrRecExtraVars; ++iRecVar)
  {
    cf_assert(nbrCellExtraVars <= recCoefs[iRecVar].size());

    // dereference extra variables to be reconstructed
    RealVector& recExtraVar = *recExtraVars[iRecVar];

    // reconstruct the extra variables
    recExtraVar = 0.0;
    for (CFuint iCellVar = 0; iCellVar < nbrCellExtraVars; ++iCellVar)
    {
      recExtraVar += recCoefs[iRecVar][iCellVar]*(*cellExtraVars[iCellVar]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFV::reconstructGradients(const vector< vector< RealVector >* >& cellGradients,
                                                       vector< vector< RealVector* > >& recGradients,
                                                       const vector< vector< CFreal > >& recCoefs,
                                                       const CFuint nbrCellGrads)
{
  // number of gradients to be reconstructed
  const CFuint nbrRecGrads = recCoefs.size();
  cf_assert(nbrRecGrads <= recGradients.size());
  cf_assert(nbrCellGrads <= cellGradients.size());

  // number of gradient variables
  cf_assert(nbrCellGrads > 0);
  const CFuint nbrGradVars = cellGradients[0]->size();

  // reconstruct the gradients
  for (CFuint iRecGrad = 0; iRecGrad < nbrRecGrads; ++iRecGrad)
  {
    cf_assert(nbrCellGrads <= recCoefs[iRecGrad].size());
    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      // dereference gradient to be reconstructed
      RealVector& recGrad = *recGradients[iRecGrad][iGradVar];

      // reconstruct the gradient
      recGrad = 0.0;
      for (CFuint iCV = 0; iCV < nbrCellGrads; ++iCV)
      {
        recGrad += recCoefs[iRecGrad][iCV]*(*cellGradients[iCV])[iGradVar];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ReconstructStatesSpectralFV::reconstructPhysVar
                                              (const CFuint iVar,
                                               const std::vector< Framework::State* >& cellStates,
                                               std::vector< Framework::State* >& recStates,
                                               const std::vector< std::vector< CFreal > >& recCoefs,
                                               const CFuint nbrCellStates)
{
  // number of reconstructed states
  const CFuint nbrRecStates = recCoefs.size();
  cf_assert(nbrRecStates <= recStates.size());
  cf_assert(nbrCellStates <= cellStates.size());

  // backup and reconstruct the physical variable
  for (CFuint iRecState = 0; iRecState < nbrRecStates; ++iRecState)
  {
    // dereference state
    State& recState = *recStates[iRecState];

    // reconstruct
    cf_assert(nbrCellStates <= recCoefs[iRecState].size());

    // reconstruct the state
    recState[iVar] = 0.0;
    for (CFuint iCellState = 0; iCellState < nbrCellStates; ++iCellState)
    {
      recState[iVar] += recCoefs[iRecState][iCellState]*(*cellStates[iCellState])[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

