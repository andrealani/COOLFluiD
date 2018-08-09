#include "Framework/MethodStrategyProvider.hh"
#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCMirrorWallHartmann.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<BCMirrorWallHartmann, FluxReconstructionSolverData, BCStateComputer, FluxReconstructionMultiFluidMHDModule >
    BCMirrorWallHartmannProvider("BCMirrorWallHartmann");

//////////////////////////////////////////////////////////////////////////////

void BCMirrorWallHartmann::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

BCMirrorWallHartmann::BCMirrorWallHartmann
(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL)
{
  CFAUTOTRACE;
}
      
//////////////////////////////////////////////////////////////////////////////

BCMirrorWallHartmann::~BCMirrorWallHartmann()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorWallHartmann::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some physical data from the model
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    
    // normal
    const RealVector& normal = normals[iState];

    CFuint nx = normal[XX];
    CFuint ny = normal[YY];
    CFuint nz = 0;
    const CFuint invFaceLength = 1./sqrt(nx*nx + ny*ny);
    nx *= invFaceLength;
    ny *= invFaceLength;

    const CFuint bn = (m_intSolPhysData)[0]*nx + (m_intSolPhysData)[1]*ny;
    const CFuint en = (m_intSolPhysData)[3]*nx + (m_intSolPhysData)[4]*ny;  
    
    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 4);
    cf_assert(ghostState.size() == 4);

    // set the physical data starting from the inner state
    m_varSet->computePhysicalData(intState,m_intSolPhysData);

  (m_ghostSolPhysData)[0] = (m_intSolPhysData)[0]; //- 2*bn*nx;	//Bx
  (m_ghostSolPhysData)[1] = (m_intSolPhysData)[1]; //- 2*bn*ny;	//By
  (m_ghostSolPhysData)[2] = (m_intSolPhysData)[2]; // - 2*bn*nz;	//Bz
  (m_ghostSolPhysData)[3] = -(m_intSolPhysData)[3] - 2*en*nx;	//Ex
  (m_ghostSolPhysData)[4] = -(m_intSolPhysData)[4] - 2*en*ny;	//Ey
  (m_ghostSolPhysData)[5] = -(m_intSolPhysData)[5] - 2*en*nz;	//Ez
  (m_ghostSolPhysData)[6] = (m_intSolPhysData)[6];			//Psi
  (m_ghostSolPhysData)[7] = -(m_intSolPhysData)[7];			//Phi
 } 
///MultiFluidMHD mirror Condition in 2D
  const CFuint endEM = 8;
 
  //set the densities
  for (CFuint i = 0 ; i < nbrStates; i++){
    (m_ghostSolPhysData)[endEM + i] = (m_intSolPhysData)[endEM + i];
  }
 
  //set the Velocities
  for (CFuint i = 0 ; i < nbrStates; i++){
    
        // normal
    const RealVector& normal = normals[i];

    CFuint nx = normal[XX];
    CFuint ny = normal[YY];
    CFuint nz = 0;
    const CFuint invFaceLength = 1./sqrt(nx*nx + ny*ny);
    nx *= invFaceLength;
    ny *= invFaceLength;
 
    CFreal un_i = (m_intSolPhysData)[endEM + nbrStates + 2*i]*nx + (m_intSolPhysData)[endEM + nbrStates + 2*i + 1]*ny; 
    (m_ghostSolPhysData)[endEM + nbrStates + 2*i] = (m_intSolPhysData)[endEM + nbrStates + 2*i] - 2*un_i*nx;
    (m_ghostSolPhysData)[endEM + nbrStates + 2*i + 1] = (m_intSolPhysData)[endEM + nbrStates + 2*i + 1] - 2*un_i*ny;
  } 
 
  //set the Temperatures
  for (CFuint i = 0 ; i < nbrStates; i++){
    (m_ghostSolPhysData)[endEM + nbrStates + 2*nbrStates + i] = (m_intSolPhysData)[endEM + nbrStates + 2*nbrStates + i];
    cf_assert((m_intSolPhysData)[endEM + nbrStates + 2*nbrStates + i] > 0.);
  } 
 
}
//////////////////////////////////////////////////////////////////////////////

void BCMirrorWallHartmann::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 5);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
//     RealVector& pGradI = *intGrads  [iState][0];
//     RealVector& pGradG = *ghostGrads[iState][0];
//     const CFreal nPGrad = pGradI[XX]*normal[XX] + pGradI[YY]*normal[YY] + pGradI[ZZ]*normal[ZZ];
//     pGradG = pGradI - 2.0*nPGrad*normal;
    *ghostGrads[iState][0] = *intGrads[iState][0];

    // velocity
    *ghostGrads[iState][1] = *intGrads[iState][1];
    *ghostGrads[iState][2] = *intGrads[iState][2];
    *ghostGrads[iState][3] = *intGrads[iState][3];
  

/*    bool m_heatflux = true;

    if (m_heatFlux)
    {
      // temperature
      RealVector& tempGradI = *intGrads  [iState][4];
      RealVector& tempGradG = *ghostGrads[iState][4];
      const CFreal nTempGrad = tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY];
      tempGradG = tempGradI - 2.0*nTempGrad*normal + m_wallQ*normal;
    }
    else
    {
      *ghostGrads[iState][4] = *intGrads[iState][4];
    }
*/

  }
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorWallHartmann::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();
  
  m_varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in BCMirrorWallHartmann!");
  }

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

