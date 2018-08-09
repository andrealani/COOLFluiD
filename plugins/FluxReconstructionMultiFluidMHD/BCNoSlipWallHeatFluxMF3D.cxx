#include "Framework/MethodStrategyProvider.hh"

#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DVarSet.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCNoSlipWallHeatFluxMF3D.hh"

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

  namespace Physics {
    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCNoSlipWallHeatFluxMF3D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMultiFluidMHDModule >
  BCNoSlipWallHeatFluxMF3DProvider("NoSlipWallHeatFluxMF3D");

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallHeatFluxMF3D::BCNoSlipWallHeatFluxMF3D(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);

  m_wallT = 0.0; //std::vector<CFreal>();
   setParameter("T",&m_wallT);

  m_wallQ = 0.0; // std::vector<CFreal>();
   setParameter("q",&m_wallQ);

  m_heatFlux= true;
   setParameter("HeatFlux",&m_heatFlux);
   
  _Ez0 = -0.1;
  setParameter("Ez0",&_Ez0);

  m_changeToIsoT = MathTools::MathConsts::CFuintMax();
   setParameter("ChangeToIsoT",&m_changeToIsoT);
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallHeatFluxMF3D::~BCNoSlipWallHeatFluxMF3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxMF3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Ez0","Imposed Ez0");
  options.addConfigOption< CFreal/*std::vector<CFreal>*/ >("T","wall static temperature");
  options.addConfigOption< CFreal/*std::vector<CFreal>*/ >("q","wall heat flux");
  options.addConfigOption< bool >("HeatFlux","bool to tell if the wall has constant heat flux, default true.");
  options.addConfigOption< CFuint >("ChangeToIsoT","Iteration after which to switch to an isothermal BC.");
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxMF3D::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  const CFuint nbSpecies = m_varSet->getModel()->getNbScalarVars(0);

  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());
  
  ///Maxwell Equations Perfectly Conducting Wall Condition
  const CFreal muZero = 1.2566370614359e-6;
  const CFreal omega = 1.82e5; //14857.14286;
  const CFreal Ez0 = _Ez0;  

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  if (iter >= m_changeToIsoT && m_heatFlux)
  {
    m_heatFlux = false;
  }

  // get some physical data from the model
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);

//   const CFreal dy = ghostState->getCoordinates()[YY] - innerState->getCoordinates()[YY];
  const CFreal dr = MathFunctions::getDistance(ghostState.getCoordinates(),
						  intState.getCoordinates());
  //cout<<"dr = "<< dr <<"\n";

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

  (m_ghostSolPhysData)[0] = (m_intSolPhysData)[0] - dr*muZero*omega*Ez0;	//Bx
  (m_ghostSolPhysData)[1] = -(m_intSolPhysData)[1]; //By
  (m_ghostSolPhysData)[2] = (m_intSolPhysData)[2]; //Bz
  (m_ghostSolPhysData)[3] = -(m_intSolPhysData)[3]; //Ex
  (m_ghostSolPhysData)[4] = -(m_intSolPhysData)[4]; //Ey
  (m_ghostSolPhysData)[5] = -(m_intSolPhysData)[5]; //Ez
  (m_ghostSolPhysData)[6] = -(m_intSolPhysData)[6];	//Psi
  (m_ghostSolPhysData)[7] = -(m_intSolPhysData)[7];	//Phi

     ///MultiFluidMHD No Slip Isothermal Condition in 2D
     // here a fix is needed in order to have always ghostT > 0
     // if ghostT < 0  then the inner value is set
     const CFuint endEM = 8;
 
    
    if (m_heatFlux)
    {
  
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
       const CFreal T = _wallTemp[i];
       (m_ghostSolPhysData)[endEM + nbSpecies + 2*nbSpecies + i] = 2.*T - (m_intSolPhysData)[endEM + nbSpecies + 2*nbSpecies + i];
       cf_assert(2.*T - (m_intSolPhysData)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
     } 

/*      // set the physical data for the ghost state
      m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::VX]  = -m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VY]  = -m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::VZ]  = -m_intSolPhysData[EulerTerm::VZ];// negate velocity )-- average = 0
      m_ghostSolPhysData[EulerTerm::V] = sqrt(m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX]+m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY]+m_ghostSolPhysData[EulerTerm::VZ]*m_ghostSolPhysData[EulerTerm::VZ]);
      m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P];
      m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
                                            + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                                  m_intSolPhysData[EulerTerm::V]*
                                                  m_intSolPhysData[EulerTerm::V]
                                         )/m_ghostSolPhysData[EulerTerm::RHO];
      m_ghostSolPhysData[EulerTerm::A] = sqrt(gamma*m_intSolPhysData[EulerTerm::P]/m_intSolPhysData[EulerTerm::RHO]);
      m_ghostSolPhysData[EulerTerm::T] = m_intSolPhysData[EulerTerm::T];
*/  
  }
    else
    {
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
       const CFreal T = _wallTemp[i];
       (m_ghostSolPhysData)[endEM + nbSpecies + 2*nbSpecies + i] = (m_intSolPhysData)[endEM + nbSpecies + 2*nbSpecies + i];
       //(*ghostState)[endEM + nbSpecies + 2*nbSpecies + i] = 2.*Twall - (m_intSolPhysData)[endEM + nbSpecies + 2*nbSpecies + i];
       cf_assert(2.*T - (m_intSolPhysData)[endEM + nbSpecies + 2*nbSpecies + i] > 0.); 
     } 

      // set the ghost state from its physical data
      m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxMF3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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
    
    if (m_heatFlux)
    {
      // temperature
      RealVector& tempGradI = *intGrads  [iState][4];
      RealVector& tempGradG = *ghostGrads[iState][4];
      const CFreal nTempGrad = tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY] + tempGradI[ZZ]*normal[ZZ];
      //tempGradG = tempGradI - 2.0*nTempGrad*normal + m_wallQ*normal;
    }
    else
    {
      *ghostGrads[iState][4] = *intGrads[iState][4];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxMF3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get Euler 2D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo< MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();

  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in BCNoSlipWallHeatFluxMF3D!");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

