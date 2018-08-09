#include "MultiFluidMHD/DiffMFMHD3DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCSubOutletPPCWRhoiViTi3D.hh"

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
    BCSubOutletPPCWRhoiViTi3D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMultiFluidMHDModule >
  BCSubOutletPPCWRhoiViTi3DProvider("BCSubOutletPPCWRhoiViTi3D");

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletPPCWRhoiViTi3D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

BCSubOutletPPCWRhoiViTi3D::BCSubOutletPPCWRhoiViTi3D(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_pressure()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_pressure = std::vector<CFreal>();       //1.0;
  setParameter("P",&m_pressure);
}

//////////////////////////////////////////////////////////////////////////////

BCSubOutletPPCWRhoiViTi3D::~BCSubOutletPPCWRhoiViTi3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletPPCWRhoiViTi3D::computeGhostStates(const vector< State* >& intStates,
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

    CFreal nx = normal[XX];
    CFreal ny = normal[YY];
    CFreal nz = 0;
    const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
    nx *= invFaceLength;
    ny *= invFaceLength;

  cf_assert(m_varSet.isNotNull());
  
  const CFreal bn = (m_intSolPhysData)[0]*nx + (m_intSolPhysData)[1]*ny + (m_intSolPhysData)[2]*nz;
  const CFreal en = (m_intSolPhysData)[3]*nx + (m_intSolPhysData)[4]*ny + (m_intSolPhysData)[5]*nz;  
//  const CFreal chi = m_varSet->getModel()->getDivECleaningConst();

  (m_ghostSolPhysData)[0] = (m_intSolPhysData)[0] - 2*bn*nx;	//Bx
  (m_ghostSolPhysData)[1] = (m_intSolPhysData)[1] - 2*bn*ny;	//By
  (m_ghostSolPhysData)[2] = (m_intSolPhysData)[2] - 2*bn*nz;	//Bz
  (m_ghostSolPhysData)[3] = -(m_intSolPhysData)[3] + 2*en*nx;	//Ex
  (m_ghostSolPhysData)[4] = -(m_intSolPhysData)[4] + 2*en*ny;	//Ey
  (m_ghostSolPhysData)[5] = -(m_intSolPhysData)[5] + 2*en*nz;	//Ez
  (m_ghostSolPhysData)[6] = (m_intSolPhysData)[6];			//Psi
  (m_ghostSolPhysData)[7] = -(m_intSolPhysData)[7];			//Phi
  
}
///MultiFluidMHD Subsonic Outlet Condition imposing pressure in 3D
 // here a fix is needed in order to have always ghostT > 0
 // if ghostT < 0  then the inner value is set
  const CFuint endEM = 8;
  const CFuint firstSpecies = m_varSet->getModel()->getFirstScalarVar(0);  
  const CFuint firstVelocity = m_varSet->getModel()->getFirstScalarVar(1);   
  const CFuint firstTemperature = m_varSet->getModel()->getFirstScalarVar(2); 

  CFuint dim = 3;
 
  (m_ghostSolPhysData)[endEM] = (m_intSolPhysData)[endEM]; 							//RHO

  
  for (CFuint ie = 0; ie < nbrStates; ++ie) {
    const CFreal P = m_pressure[ie];

    // dereference states
    State& intState   = (*intStates[ie]);
    State& ghostState = (*ghostStates[ie]);

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

    // set the physical data starting from the inner state
    m_varSet->computePhysicalData(intState,m_intSolPhysData);


    (m_ghostSolPhysData)[firstSpecies + ie] = (m_intSolPhysData)[firstSpecies + ie];			//yi
    (m_ghostSolPhysData)[firstVelocity + dim*ie] = (m_intSolPhysData)[firstVelocity + dim*ie];			//Vxi
    (m_ghostSolPhysData)[firstVelocity + dim*ie + 1] = (m_intSolPhysData)[firstVelocity + dim*ie + 1];		//Vyi
    (m_ghostSolPhysData)[firstVelocity + dim*ie + 2] = (m_intSolPhysData)[firstVelocity + dim*ie + 2];		//Vzi
    (m_ghostSolPhysData)[firstTemperature + 4*ie] = (m_intSolPhysData)[firstTemperature + 4*ie];		//Ti
    (m_ghostSolPhysData)[firstTemperature + 4*ie + 1] = 2*P - 
						    (m_intSolPhysData)[firstTemperature + 4*ie + 1];	//Pi
						    
    const CFreal Vi2 = (m_ghostSolPhysData)[firstVelocity + dim*ie]*(m_ghostSolPhysData)[firstVelocity + dim*ie] +
			(m_ghostSolPhysData)[firstVelocity + dim*ie + 1]*(m_ghostSolPhysData)[firstVelocity + dim*ie + 1] + 
            (m_ghostSolPhysData)[firstVelocity + dim*ie + 2]*(m_ghostSolPhysData)[firstVelocity + dim*ie + 2];
    const CFreal rhoi =(m_ghostSolPhysData)[endEM]*(m_ghostSolPhysData)[firstSpecies + ie];
    
    (m_ghostSolPhysData)[firstTemperature + 4*ie + 2] = sqrt(gamma*(m_ghostSolPhysData)[firstTemperature + 4*ie + 1]/rhoi);		//ai
    (m_ghostSolPhysData)[firstTemperature + 4*ie + 3] = (0.5*rhoi*Vi2 + 
						    gammaDivGammaMinus1*(m_ghostSolPhysData)[firstTemperature + 4*ie + 1])/rhoi;	//Hi


    m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);

  }
    
}

//////////////////////////////////////////////////////////////////////////////

void BCSubOutletPPCWRhoiViTi3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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

void BCSubOutletPPCWRhoiViTi3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  m_varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> >();
  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData);
  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  
//   const CFuint nbSpecies = _updateVarSet->getModel()->getNbScalarVars(0);  
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     cf_assert(_wallTemp[i] >= 0.0);
//   }

  // adimensionalize the temperature
//   _wallTemp /= _updateVarSet->getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

