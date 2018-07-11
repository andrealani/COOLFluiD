#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BCSuperOutlet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSuperOutlet,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionModule >
  BCSuperOutletProvider("SuperOutlet");

//////////////////////////////////////////////////////////////////////////////

BCSuperOutlet::BCSuperOutlet(const std::string& name) :
  BCStateComputer(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCSuperOutlet::~BCSuperOutlet()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::computeGhostStates(const vector< State* >& intStates,
                                       vector< State* >& ghostStates,
                                       const std::vector< RealVector >& normals,
                                       const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates == ghostStates.size());
  cf_assert(nbrStates == normals.size());

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    *ghostStates[iState] = *intStates[iState];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                          std::vector< std::vector< RealVector* > >& ghostGrads,
                                          const std::vector< RealVector >& normals,
                                          const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);
  const CFuint nbrGradVars = intGrads[0].size();

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    // normal
    const RealVector& normal = normals[iState];

    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      const RealVector& varGradI =  *intGrads[iState][iGradVar];
      RealVector& varGradG =  *ghostGrads[iState][iGradVar];
      const CFreal nVarGrad = MathTools::MathFunctions::innerProd(varGradI, normal);
      varGradG = varGradI - 2.0*nVarGrad*normal;

//      *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar]; //0;//
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::unsetup()
{
  CFAUTOTRACE;

  // unsetup of the parent class
  BCStateComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

