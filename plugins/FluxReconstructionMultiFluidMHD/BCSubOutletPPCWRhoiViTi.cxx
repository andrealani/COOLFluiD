#include "Framework/MethodStrategyProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/MeshData.hh"
#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCSubOutletPPCWRhoiViTi.hh"

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
    BCSubOutletPPCWRhoiViTi,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMultiFluidMHDModule >
  BCSubOutletPPCWRhoiViTiProvider("BCSubOutletPPCWRhoiViTi");

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletPPCWRhoiViTi::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFreal> >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

BCSubOutletPPCWRhoiViTi::BCSubOutletPPCWRhoiViTi(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_pressure()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_pressure = std::vector<CFreal>(); //1.0;
  setParameter("P",&m_pressure);
}

//////////////////////////////////////////////////////////////////////////////

BCSubOutletPPCWRhoiViTi::~BCSubOutletPPCWRhoiViTi()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletPPCWRhoiViTi::computeGhostStates(const vector< State* >& intStates,
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

//  CFLog(VERBOSE, "\n\n\n Error 1 \n");
 
///MultiFluidMHD Subsonic Outlet Condition imposing pressure in 2D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
  const CFuint endEM = 8;
  const CFuint firstSpecies = m_varSet->getModel()->getFirstScalarVar(0);  
  const CFuint firstVelocity = m_varSet->getModel()->getFirstScalarVar(1);   
  const CFuint firstTemperature = m_varSet->getModel()->getFirstScalarVar(2); 
/*
    CFLog(VERBOSE, "\n\n\n\n\n m_intSolPhysData size  = " << m_intSolPhysData.size() << "\n\n\n\n\n"); 
    CFLog(VERBOSE, "\n\n\n\n\n First velocity = " << firstVelocity << "\n\n\n\n\n");
    CFLog(VERBOSE, "\n\n\n\n\n First temperature = " << firstTemperature << "\n\n\n\n\n");
    CFLog(VERBOSE, "\n\n\n\n\n nbspecies = " << nbSpecies << "\n\n\n\n\n");
*/
  //const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal Boltz = m_varSet->getModel()->getK();
  const CFreal molecularMass1 = m_varSet->getModel()->getMolecularMass1();
  const CFreal molecularMass2 = m_varSet->getModel()->getMolecularMass2();
  const CFreal molecularMass3 = m_varSet->getModel()->getMolecularMass3();
  const CFreal R1 = Boltz/molecularMass1;
  const CFreal R2 = Boltz/molecularMass2;
  const CFreal R3 = Boltz/molecularMass3;

  (m_ghostSolPhysData)[endEM] = (m_intSolPhysData)[endEM]; 							//RHO

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  { 

      // dereference states
      State& intState   = (*intStates[iState]);
      State& ghostState = (*ghostStates[iState]);

      // set the physical data starting from the inner state
      m_varSet->computePhysicalData(intState, m_intSolPhysData);      

    // normal
    const RealVector& normal = normals[iState];

    CFreal nx = normal[XX];
    CFreal ny = normal[YY];
    CFreal nz = 0;
    const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
    nx *= invFaceLength;
    ny *= invFaceLength;

    cf_assert(m_varSet.isNotNull());

    const CFreal bn = (m_intSolPhysData)[0]*nx + (m_intSolPhysData)[1]*ny;
    const CFreal en = (m_intSolPhysData)[3]*nx + (m_intSolPhysData)[4]*ny;  


 // CFLog(VERBOSE, "\n\n\n Error 2 \n");

  (m_ghostSolPhysData)[0] = (m_intSolPhysData)[0] /*+ 2*bn*nx*/;	//Bx
  (m_ghostSolPhysData)[1] = (m_intSolPhysData)[1] /*+ 2*bn*ny*/;	//By
  (m_ghostSolPhysData)[2] = (m_intSolPhysData)[2] /*+ 2*bn*nz*/;	//Bz
  (m_ghostSolPhysData)[3] = (m_intSolPhysData)[3] /*- 2*en*nx*/;	//Ex
  (m_ghostSolPhysData)[4] = (m_intSolPhysData)[4] /*- 2*en*ny*/;	//Ey
  (m_ghostSolPhysData)[5] = (m_intSolPhysData)[5] /*- 2*en*nz*/;	//Ez
  (m_ghostSolPhysData)[6] = (m_intSolPhysData)[6];			//Psi
  (m_ghostSolPhysData)[7] = (m_intSolPhysData)[7];			//Phi

//    cf_assert(intState.size() == 4);
//    cf_assert(ghostState.size() == 4);

/*
  CFLog(VERBOSE, "\n\n\n Error 3 \n");

    CFLog(VERBOSE, "\n\n\n\n\n nbSpecies Value = " << nbSpecies << "\n\n\n\n\n");
    CFLog(VERBOSE, "\n\n\n\n\n m_intSolPhysData size Value = " << m_intSolPhysData << "\n\n\n\n\n");
*/

  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
   const CFreal P = m_pressure[ie];
    (m_ghostSolPhysData)[firstSpecies + ie] = (m_intSolPhysData)[firstSpecies          + ie];			//yi                 
    (m_ghostSolPhysData)[firstVelocity + 2*ie] = (m_intSolPhysData)[firstVelocity + 2*ie];			//Vxi
    (m_ghostSolPhysData)[firstVelocity + 2*ie + 1] = (m_intSolPhysData)[firstVelocity + 2*ie + 1];		//Vyi
    (m_ghostSolPhysData)[firstTemperature + 4*ie] = (m_intSolPhysData)[firstTemperature + 4*ie];		//Ti
    (m_ghostSolPhysData)[firstTemperature + 4*ie + 1] = 2*P - 
						    (m_intSolPhysData)[firstTemperature + 4*ie + 1];	//P
						    
    const CFreal Vi2 = (m_ghostSolPhysData)[firstVelocity + 2*ie]*(m_ghostSolPhysData)[firstVelocity + 2*ie] +
			(m_ghostSolPhysData)[firstVelocity + 2*ie + 1]*(m_ghostSolPhysData)[firstVelocity + 2*ie + 1];



    const CFreal rhoi =(m_ghostSolPhysData)[endEM]*(m_ghostSolPhysData)[firstSpecies + ie];

    (m_ghostSolPhysData)[firstTemperature + 4*ie + 2] = sqrt(gamma*(m_ghostSolPhysData)[firstTemperature + 4*ie + 1]/rhoi);		//ai
    (m_ghostSolPhysData)[firstTemperature + 4*ie + 3] = (0.5*rhoi*Vi2 + 
						    gammaDivGammaMinus1*(m_ghostSolPhysData)[firstTemperature + 4*ie + 1])/rhoi;	//Hi

     // Temporary fix to one fluid 
  (m_ghostSolPhysData)[endEM] = (m_ghostSolPhysData)[firstTemperature + 1]/(R1*(m_ghostSolPhysData)[firstTemperature]);                //RHO is computed with the pressure

  (m_ghostSolPhysData)[firstSpecies] = 1;											      //y = 1
  const CFreal V2 = (m_ghostSolPhysData)[firstVelocity]*(m_ghostSolPhysData)[firstVelocity] +
                      (m_ghostSolPhysData)[firstVelocity + 1]*(m_ghostSolPhysData)[firstVelocity + 1];			      //Vx and Vy are correct;
  const CFreal rho =(m_ghostSolPhysData)[endEM];										      // we use this rho
  (m_ghostSolPhysData)[firstTemperature + 2] = sqrt(gamma*(m_ghostSolPhysData)[firstTemperature + 1]/rho);                         //ai
  (m_ghostSolPhysData)[firstTemperature + 3] = (0.5*rho*V2 +
                                            gammaDivGammaMinus1*(m_ghostSolPhysData)[firstTemperature + 1])/rho;                //Hi


  }
   m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState); 

   }      
}

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletPPCWRhoiViTi::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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
}

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletPPCWRhoiViTi::setup()
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
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in BCSuboutletPPCWRhoiViTi!");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData);
  m_varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
 
  // non-dimensionalize pressure and temperature
  //m_pressure /= m_varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

