#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/BCSuperInletProjMHD.hh"

#include "Common/NotImplementedException.hh"

#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSuperInletProjMHD,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMHDModule >
  BCSuperInletProjMHDProvider("SuperInletProjMHD");

//////////////////////////////////////////////////////////////////////////////

BCSuperInletProjMHD::BCSuperInletProjMHD(const std::string& name) :
  BCStateComputer(name),
  m_initialSolutionMap()
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);

  m_rhoBC = 1.0;
  setParameter("RhoBC",&m_rhoBC);
  
  m_pBC = 0.108;
  setParameter("pBC",&m_pBC);
}

//////////////////////////////////////////////////////////////////////////////

BCSuperInletProjMHD::~BCSuperInletProjMHD()
{
  CFAUTOTRACE;
  
  if (m_initialSolutionMap.size() > 0) 
  {
    for (CFuint i = 0; i < m_initialSolutionMap.size(); ++i) 
    {
      deletePtr(m_initialSolutionMap[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletProjMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("RhoBC","Boundary rho value.");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("pBC","Boundary p value.");
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletProjMHD::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());
  
  bool is3D = (normals[0]).size()==3;

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 9 || 7);
    cf_assert(ghostState.size() == 9 || 7);

    // during initialization phase we store the initial solution values to be used as BC
    // This BC (the PFSS solution fixed on the inlet) was suggested by Jon Linker and
    // is evaluated by default. When using Jens' or Dana's BCs for the magnetic field
    // the B-field values in the ghost cells are overwritten.

    // initialize PFSS B solution vector
    RealVector B_PFSS_dimless(3);

    // get PFSS B solution from m_initialSolutionMap
    if (m_initialSolutionMap.size() > 0) 
    {
      /// map faces to corresponding TRS and index inside that TRS
      SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
      const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(m_face->getID());
      const string name = m_trsNames[0];//getCurrentTRS()->getName();
      SafePtr<RealVector> initialValues = m_initialSolutionMap.find(name);
      const CFuint nbVars = is3D ? 3 : 2;//m_initialSolutionIDs.size();
      const CFuint startID = faceIdx*nbVars;
      for (CFuint i = 0; i < nbVars; ++i) 
      {
        const CFuint varID = is3D ? i+4 : i+3;//m_initialSolutionIDs[i];
        const CFuint idx = startID+i;
        cf_assert(idx < initialValues->size());

        // save the Poisson PFSS solution in a vector
        B_PFSS_dimless[i] = (*initialValues)[idx]; 

        //(*ghostState)[varID] = 2.*(*initialValues)[idx] - (*innerState)[varID];
        //CFLog(DEBUG_MIN, "SuperInletProjectionParallel5::setGhostState() => [" << varID << "] => " << (*initialValues)[idx] << " | " << (*innerState)[varID] << "\n");
      }
    }
    
    const RealVector& normal = normals[iState];
      
    const CFreal xI_dimless = coords[iState][XX];
    const CFreal yI_dimless = coords[iState][YY];
    const CFreal zI_dimless = coords[iState][ZZ];
    const CFreal rI_dimless = (coords[iState]).norm2();
    
    // to be checked if imposed minimum is needed (suspected problems at poles otherwise)
    const CFreal rhoI_dimless = std::max(std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless),1.0e-4);
    
    // density
    ghostState[0] = 2.0*m_rhoBC - intState[0];
    
    // p
    ghostState[7] = 2.0*m_pBC - intState[7];

    // phi
    ghostState[8] = intState[8];
    
    // V
    ghostState[1] = -intState[1];
    ghostState[2] = -intState[2];
    ghostState[3] = -intState[3]; 
    
    // B
    // Br bnd from PFSS solution
    const CFreal BrBoundary_dimless = xI_dimless/rhoI_dimless*B_PFSS_dimless[0] + yI_dimless/rI_dimless*B_PFSS_dimless[1] + zI_dimless/rI_dimless*B_PFSS_dimless[2];

    const CFreal BxI_dimless = intState[4];
    const CFreal ByI_dimless = intState[5];
    const CFreal BzI_dimless = intState[6];
    const CFreal BrI_dimless = xI_dimless/rI_dimless*BxI_dimless + yI_dimless/rI_dimless*ByI_dimless + zI_dimless/rI_dimless*BzI_dimless;
    const CFreal BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*ByI_dimless - rhoI_dimless/rI_dimless*BzI_dimless;
    const CFreal BphiI_dimless = -yI_dimless/rhoI_dimless*BxI_dimless + xI_dimless/rhoI_dimless*ByI_dimless;

    // B bnd, Btheta and Bphi bnd are set equal to inner values
    const CFreal BxBoundary_dimless = xI_dimless/rI_dimless*BrBoundary_dimless - yI_dimless/rhoI_dimless*BphiI_dimless + xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BthetaI_dimless;
    const CFreal ByBoundary_dimless = yI_dimless/rI_dimless*BrBoundary_dimless + xI_dimless/rhoI_dimless*BphiI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BthetaI_dimless;
    const CFreal BzBoundary_dimless = zI_dimless/rI_dimless*BrBoundary_dimless - rhoI_dimless/rI_dimless*BthetaI_dimless;

    ghostState[4] = 2.0*BxBoundary_dimless - BxI_dimless;
    ghostState[5] = 2.0*ByBoundary_dimless - ByI_dimless;
    ghostState[6] = 2.0*BzBoundary_dimless - BzI_dimless;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletProjMHD::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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

void BCSuperInletProjMHD::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = true;

//  // get MHD 3D varset
//  m_varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
//  if (m_varSet.isNull())
//  {
//    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHD3DProjectionVarSet in BCSuperInletProjMHD!");
//  }
//
//  // resize the physical data for internal and ghost solution points
//  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
//  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

