#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionTurb/DiffBndCorrectionsRHSJacobFluxReconstructionTurb.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionTurb, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionNavierStokesModule >
DiffBndCorrectionsRHSJacobTurbFluxReconstructionProvider("DiffBndCorrectionsRHSJacobTurb");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionTurb::DiffBndCorrectionsRHSJacobFluxReconstructionTurb(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstructionNS(name),
  socket_wallDistance("wallDistance"),
  m_closestSolToFlxIdx(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionTurb::~DiffBndCorrectionsRHSJacobFluxReconstructionTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionTurb::setup()
{
  DiffBndCorrectionsRHSJacobFluxReconstructionNS::setup();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  // get closest sol indices
  m_closestSolToFlxIdx = frLocalData[0]->getClosestSolToFlxIdx();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionTurb::unsetup()
{
  DiffBndCorrectionsRHSJacobFluxReconstructionNS::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    DiffBndCorrectionsRHSJacobFluxReconstructionTurb::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffBndCorrectionsRHSJacobFluxReconstructionNS::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionTurb::computeInterfaceFlxCorrection()
{ 
  // Get the wall distance
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  // compute the riemann flux in the flx pnts
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // local flux point indices in the left and right cell
    const CFuint flxPntIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
    const CFuint closestSolIdx = (*m_closestSolToFlxIdx)[flxPntIdx];
    
    const CFuint stateID = (*m_cellStates)[closestSolIdx]->getLocalID();
    
    // Set the wall distance before computing the turbulent viscosity
    navierStokesVarSet->setWallDistance(wallDist[stateID]);
      
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_flxPntGhostGrads[iFlxPnt][iVar]) + *(m_cellGradFlxPnt[iFlxPnt][iVar]))/2.0;

      m_avgSol[iVar] = ((*(m_flxPntGhostSol[iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[iFlxPnt]))[iVar])/2.0; 
    }
    // prepare the flux computation
    prepareFluxComputation();

    // compute FI
    //m_flxPntRiemannFlux[iFlxPnt] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0);
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
  
    // compute FI in the local frame
    m_cellFlx[iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
