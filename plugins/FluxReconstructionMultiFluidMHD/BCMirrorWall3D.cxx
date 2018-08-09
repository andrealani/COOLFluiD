#include "Framework/MethodStrategyProvider.hh"

#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCMirrorWall3D.hh"

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

Framework::MethodStrategyProvider<
    BCMirrorWall3D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMultiFluidMHDModule >
  BCMirrorWall3DProvider("BCMirrorWall3D");

//////////////////////////////////////////////////////////////////////////////

void BCMirrorWall3D::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

BCMirrorWall3D::BCMirrorWall3D(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCMirrorWall3D::~BCMirrorWall3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorWall3D::computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
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
    const CFuint nbSpecies = m_varSet->getModel()->getNbScalarVars(0);

    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

    // set the physical data starting from the inner state
    m_varSet->computePhysicalData(intState,m_intSolPhysData);

    // normal
    const RealVector& normal = normals[iState];

    CFreal nx = normal[XX];
    CFreal ny = normal[YY];
    CFreal nz = normal[ZZ];
    const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
    nx *= invFaceLength;
    ny *= invFaceLength;
    nz *= invFaceLength;

  cf_assert(m_varSet.isNotNull());
  
  const CFreal bn = (m_intSolPhysData)[0]*nx + (m_intSolPhysData)[1]*ny + (m_intSolPhysData)[2]*nz;
  const CFreal en = (m_intSolPhysData)[3]*nx + (m_intSolPhysData)[4]*ny + (m_intSolPhysData)[5]*nz;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (m_ghostSolPhysData)[0] = -(m_intSolPhysData)[0] + 2*bn*nx;	//Bx
  (m_ghostSolPhysData)[1] = -(m_intSolPhysData)[1] + 2*bn*ny;	//By
  (m_ghostSolPhysData)[2] = -(m_intSolPhysData)[2] + 2*bn*nz;	//Bz
  (m_ghostSolPhysData)[3] = (m_intSolPhysData)[3] - 2*en*nx;	//Ex
  (m_ghostSolPhysData)[4] = (m_intSolPhysData)[4] - 2*en*ny;	//Ey
  (m_ghostSolPhysData)[5] = (m_intSolPhysData)[5] - 2*en*nz;	//Ez
  (m_ghostSolPhysData)[6] = (m_intSolPhysData)[6];			//Psi
  (m_ghostSolPhysData)[7] = -(m_intSolPhysData)[7];			//Phi

  
///MultiFluidMHD mirror Condition in 3D
  const CFuint endEM = 8;
 
  //set the densities
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (m_ghostSolPhysData)[endEM + i] = (m_intSolPhysData)[endEM + i];
  }
 
  //set the Velocities
  CFuint dim = 3;

  for (CFuint i = 0 ; i < nbSpecies; i++){

    CFreal un_i = (m_intSolPhysData)[endEM + nbSpecies + dim*i]*nx + (m_intSolPhysData)[endEM + nbSpecies + dim*i + 1]*ny + (m_intSolPhysData)[endEM + nbSpecies + dim*i + 2]*nz; 
    (m_ghostSolPhysData)[endEM + nbSpecies + dim*i]     = (m_intSolPhysData)[endEM + nbSpecies + dim*i] - 2*un_i*nx;
    (m_ghostSolPhysData)[endEM + nbSpecies + dim*i + 1] = (m_intSolPhysData)[endEM + nbSpecies + dim*i + 1] - 2*un_i*ny;
    (m_ghostSolPhysData)[endEM + nbSpecies + dim*i + 2] = (m_intSolPhysData)[endEM + nbSpecies + dim*i + 2] - 2*un_i*nz;
  } 
 
  //set the Temperatures
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (m_ghostSolPhysData)[endEM + nbSpecies + dim*nbSpecies + i] = (m_intSolPhysData)[endEM + nbSpecies + dim*nbSpecies + i];
    cf_assert((m_intSolPhysData)[endEM + nbSpecies + dim*nbSpecies + i] > 0.); 
  } 

    // set the ghost state from its physical data
    m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}
//////////////////////////////////////////////////////////////////////////////

void BCMirrorWall3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 5);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
    RealVector& presGradI = *intGrads  [iState][0];
    RealVector& presGradG = *ghostGrads[iState][0];
    const CFreal nPresGrad = (presGradI[XX]*normal[XX] + presGradI[YY]*normal[YY] + presGradI[ZZ]*normal[ZZ]);
    presGradG = presGradI - 2.0*nPresGrad*normal;

    // temperature
    RealVector& tempGradI = *intGrads  [iState][4];
    RealVector& tempGradG = *ghostGrads[iState][4];
    const CFreal nTempGrad = (tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY] + tempGradI[ZZ]*normal[ZZ]);
    tempGradG = tempGradI - 2.0*nTempGrad*normal;

    // velocity
    RealVector& uGradI = *intGrads  [iState][1];
    RealVector& uGradG = *ghostGrads[iState][1];
    RealVector& vGradI = *intGrads  [iState][2];
    RealVector& vGradG = *ghostGrads[iState][2];
    RealVector& wGradI = *intGrads  [iState][3];
    RealVector& wGradG = *ghostGrads[iState][3];

    // internal normal and tangential component
    const RealVector velocNGradI = uGradI*normal[XX] + vGradI*normal[YY] + wGradI*normal[ZZ];
    const RealVector uTGradI = uGradI - velocNGradI*normal[XX];
    const RealVector vTGradI = vGradI - velocNGradI*normal[YY];
    const RealVector wTGradI = wGradI - velocNGradI*normal[ZZ];

    // ghost normal and tangential component
    const RealVector velocNGradG = velocNGradI;
    RealVector velocTGradNI(3);
    velocTGradNI[XX] = uTGradI[XX]*normal[XX] + uTGradI[YY]*normal[YY] + uTGradI[ZZ]*normal[ZZ];
    velocTGradNI[YY] = vTGradI[XX]*normal[XX] + vTGradI[YY]*normal[YY] + vTGradI[ZZ]*normal[ZZ];
    velocTGradNI[ZZ] = wTGradI[XX]*normal[XX] + wTGradI[YY]*normal[YY] + wTGradI[ZZ]*normal[ZZ];
    const RealVector uTGradG = uTGradI - 2.0*velocTGradNI[XX]*normal;
    const RealVector vTGradG = vTGradI - 2.0*velocTGradNI[YY]*normal;
    const RealVector wTGradG = wTGradI - 2.0*velocTGradNI[ZZ]*normal;

    // compute ghost velocity gradients
    uGradG = uTGradG + velocNGradG*normal[XX];
    vGradG = vTGradG + velocNGradG*normal[YY];
    wGradG = wTGradG + velocNGradG*normal[ZZ];
  }

}


//////////////////////////////////////////////////////////////////////////////

void BCMirrorWall3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  m_varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >();
  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiFluidMHDVarSet in BCMirrorWall3D!");
  }

}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

