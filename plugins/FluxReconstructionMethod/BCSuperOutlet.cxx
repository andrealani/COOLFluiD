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
  
  addConfigOptionsTo(this);

  m_zeroGrad = false;
  setParameter("ZeroGrad",&m_zeroGrad);
}

//////////////////////////////////////////////////////////////////////////////

BCSuperOutlet::~BCSuperOutlet()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("ZeroGrad","Boolean telling whether the normal gradients should be put to zero (default false).");
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::configure ( Config::ConfigArgs& args )
{
  BCStateComputer::configure(args);
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
      if (m_zeroGrad)
      {
        const RealVector& varGradI =  *intGrads[iState][iGradVar];
        RealVector& varGradG =  *ghostGrads[iState][iGradVar];
        const CFreal nVarGrad = MathTools::MathFunctions::innerProd(varGradI, normal);
        varGradG = varGradI - 2.0*nVarGrad*normal;
      }
      else
      {
        *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar]; 
      }

//      *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar]; //0;//
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                       const std::vector< Framework::State* >& intStates,
                                       const std::vector< Framework::State* >& ghostStates,
                                       const std::vector< RealVector >& unitNormals,
                                       const std::vector< RealVector >& flxPntCoords,
                                       std::vector< RealVector* >& bndGradVars)
{
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates <= gradVarsFlxPnt.size());
  cf_assert(nbrStates <= bndGradVars.size());

  // g_b = a
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    *bndGradVars[iState] = *gradVarsFlxPnt[iState];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::computeBndStates(const std::vector< Framework::State* >& intStates,
                                     const std::vector< Framework::State* >& ghostStates,
                                     const std::vector< RealVector >& unitNormals,
                                     const std::vector< RealVector >& flxPntCoords,
                                     std::vector< RealVector* >& bndStates)
{
  const CFuint nbrStates = intStates.size();

  // U_b = U
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    *bndStates[iState] = *intStates[iState];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutlet::computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                                    std::vector< std::vector< RealVector* > >& bndGrads,
                                    const std::vector< RealVector* >& bndStates,
                                    const std::vector< RealVector >& unitNormals,
                                    const std::vector< RealVector >& flxPntCoords)
{
  // q_b = q
  copyGradients(intGrads,bndGrads);

  // zero normal gradients
  if (m_zeroGrad)
  {
    const CFuint nbrFlxPnts = intGrads.size();

    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      const CFuint nbrGradVars = intGrads[iFlx].size();

      for (CFuint iVar = 0; iVar < nbrGradVars; ++iVar)
      {
        removeNormalComponent(*bndGrads[iFlx][iVar],unitNormals[iFlx]);
      }
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

