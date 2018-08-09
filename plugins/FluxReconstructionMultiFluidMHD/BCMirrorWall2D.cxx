#include "Framework/MethodStrategyProvider.hh"

#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCMirrorWall2D.hh"
  
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
    BCMirrorWall2D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMultiFluidMHDModule >
  BCMirrorWall2DProvider("BCMirrorWall2D");

//////////////////////////////////////////////////////////////////////////////

BCMirrorWall2D::BCMirrorWall2D(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCMirrorWall2D::~BCMirrorWall2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorWall2D::computeGhostStates(const vector< State* >& intStates,
                                         vector< State* >& ghostStates,
                                         const std::vector< RealVector >& normals,
                                         const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  const CFuint nbSpecies = m_varSet->getModel()->getNbScalarVars(0);
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some physical data from the model
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {

    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

//    cf_assert(intState.size() == 4);
//    cf_assert(ghostState.size() == 4);

    // set the physical data starting from the inner state
//    m_varSet->computePhysicalData(intState,m_intSolPhysData);

    // normal
    const RealVector& normal = normals[iState];

    CFreal nx = normal[XX];
    CFreal ny = normal[YY];
    CFreal nz = 0;
    const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
    nx *= invFaceLength;
    ny *= invFaceLength;

  cf_assert(m_varSet.isNotNull());
  
  const CFreal bn = (intState)[0]*nx + (intState)[1]*ny;
  const CFreal en = (intState)[3]*nx + (intState)[4]*ny;  
//  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

  (ghostState)[0] = -(intState)[0] + 2*bn*nx;	//Bx
  (ghostState)[1] = -(intState)[1] + 2*bn*ny;	//By
  (ghostState)[2] = -(intState)[2] + 2*bn*nz;	//Bz
  (ghostState)[3] = (intState)[3] - 2*en*nx;	//Ex
  (ghostState)[4] = (intState)[4] - 2*en*ny;	//Ey
  (ghostState)[5] = (intState)[5] - 2*en*nz;	//Ez
  (ghostState)[6] = (intState)[6];			//Psi
  (ghostState)[7] = -(intState)[7];			//Phi



  //CFLog(VERBOSE, "\n\n\n\n\n Ghost State = "<< ghostState << "\n\n\n\n\n"); //got some finite values


///MultiFluidMHD mirror Condition in 2D
  const CFuint endEM = 8;
 
  //set the densities
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (ghostState)[endEM + i] = (intState)[endEM + i];
  //CFLog(VERBOSE, "\n\n\n\n\n Ghost State Densites = "<< ghostState[endEM + i] << "\n\n\n\n\n");

  }
 
  //set the Velocities
  for (CFuint i = 0 ; i < nbSpecies; i++){

    CFreal un_i = (intState)[endEM + nbSpecies + 2*i]*nx + (intState)[endEM + nbSpecies + 2*i + 1]*ny; 
    (ghostState)[endEM + nbSpecies + 2*i] = (intState)[endEM + nbSpecies + 2*i] - 2*un_i*nx;
    (ghostState)[endEM + nbSpecies + 2*i + 1] = (intState)[endEM + nbSpecies + 2*i + 1] - 2*un_i*ny;

    //CFLog(VERBOSE, "\n\n\n\n\n Ghost State Velocities = "<< ghostState[endEM + nbSpecies + 2*i] << "\n");
    //CFLog(VERBOSE, " Ghost State Velocities = "<< ghostState[endEM + nbSpecies + 2*i + 1] << "\n\n");
  } 
 
  //set the Temperatures
  for (CFuint i = 0 ; i < nbSpecies; i++){
    (ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = (intState)[endEM + nbSpecies + 2*nbSpecies + i];
    cf_assert((intState)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 

    //CFLog(VERBOSE, "\n\n\n\n\n Ghost State Temperature = "<< ghostState[endEM + nbSpecies + 2*nbSpecies + i] << "\n\n\n\n\n");

  } 

  // set the ghost state from its physical data
  //m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);

 }
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorWall2D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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
    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar];
    }
  }

/*    
    cf_assert(intGrads[iState].size() == 4);

    // normal
    const RealVector& normal = normals[iState];

    // tangential unit vector
    RealVector tangent(2);
    tangent[XX] = -normal[YY];
    tangent[YY] =  normal[XX];

    // pressure
    RealVector& presGradI = *intGrads  [iState][0];
    RealVector& presGradG = *ghostGrads[iState][0];
    const CFreal nPresGrad = (presGradI[XX]*normal[XX] + presGradI[YY]*normal[YY]);
    presGradG = presGradI - 2.0*nPresGrad*normal;

    // temperature
    RealVector& tempGradI = *intGrads  [iState][3];
    RealVector& tempGradG = *ghostGrads[iState][3];
    const CFreal nTempGrad = (tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY]);
    tempGradG = tempGradI - 2.0*nTempGrad*normal;

    // velocity
    RealVector& uGradI = *intGrads  [iState][1];
    RealVector& uGradG = *ghostGrads[iState][1];
    RealVector& vGradI = *intGrads  [iState][2];
    RealVector& vGradG = *ghostGrads[iState][2];

    // internal normal and tangential component
    const RealVector uNGradI = uGradI*normal [XX] + vGradI*normal [YY];
    const RealVector uTGradI = uGradI*tangent[XX] + vGradI*tangent[YY];

    // ghost normal and tangential component
    const RealVector uNGradG = uNGradI;
    const CFreal nGradUT = uTGradI[XX]*normal[XX] + uTGradI[YY]*normal[YY];
    const RealVector uTGradG = uTGradI - 2.0*nGradUT*normal;

    // project onto x- and y-axis
    uGradG = uNGradG*normal[XX] + uTGradG*tangent[XX];
    vGradG = uNGradG*normal[YY] + uTGradG*tangent[YY];
  }
*/
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorWall2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get Euler 2D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiFluidMHD2DVarSet in BCMirrorWall2D!");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

