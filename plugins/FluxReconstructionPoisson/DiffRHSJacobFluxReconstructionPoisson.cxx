// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionPoisson/DiffRHSJacobFluxReconstructionPoisson.hh"
#include "FluxReconstructionPoisson/FluxReconstructionPoisson.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "Poisson/PoissonDiffVarSet.hh"
#include "Poisson/PoissonConvVarSet.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Poisson;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffRHSJacobFluxReconstructionPoisson,
		       FluxReconstructionSolverData,
		       FluxReconstructionPoissonModule >
diffRHSJacobPoissonFluxReconstructionProvider("DiffRHSJacobPoisson");
  
//////////////////////////////////////////////////////////////////////////////
  
DiffRHSJacobFluxReconstructionPoisson::DiffRHSJacobFluxReconstructionPoisson(const std::string& name) :
  DiffRHSJacobFluxReconstruction(name),
  m_tempGradTermL(),
  m_tempGradTermR(),
  m_diffVarSetPoisson(CFNULL),
  m_tempStatesL(),
  m_tempStatesR(),
  socket_Bx("Bx"),
  socket_By("By"),
  socket_Bz("Bz"),
  socket_Br("Br"),
  socket_Btheta("Btheta"),
  socket_Bphi("Bphi")
{
  //addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  DiffRHSJacobFluxReconstructionPoisson::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = DiffRHSJacobFluxReconstruction::providesSockets();

  result.push_back(&socket_Bx);
  result.push_back(&socket_By);
  result.push_back(&socket_Bz);
  result.push_back(&socket_Br);
  result.push_back(&socket_Btheta);
  result.push_back(&socket_Bphi);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 33.0;
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    //for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::prepareFluxComputation()
{
//  const bool isPerturb = this->getMethodData().isPerturb();
//  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
//  m_diffVarSetPoisson->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::computeUnpertCellDiffResiduals(const CFuint side)
{
  // diffusive residual of the cell
  DiffRHSJacobFluxReconstruction::computeUnpertCellDiffResiduals(side);

  // magnetic field of the cell
  computeMagneticField(side);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::computeMagneticField(const CFuint side)
{
  DataHandle< CFreal > Bx = socket_Bx.getDataHandle();
  DataHandle< CFreal > By = socket_By.getDataHandle();
  DataHandle< CFreal > Bz = socket_Bz.getDataHandle();
  DataHandle< CFreal > Br = socket_Br.getDataHandle();
  DataHandle< CFreal > Btheta = socket_Btheta.getDataHandle();
  DataHandle< CFreal > Bphi = socket_Bphi.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const State& state = *(*m_states[side])[iSol];
    const CFuint solID = state.getLocalID();
    const RealVector& coords = state.getCoordinates();

    // B is the corrected gradient of the potential
    const RealVector& gradPhi = (*m_cellGrads[side][iSol])[0];
    const CFreal BxSolPnt = gradPhi[0];
    const CFreal BySolPnt = gradPhi[1];
    const CFreal BzSolPnt = (m_dim == 3) ? gradPhi[2] : 0.0;

    Bx[solID] = BxSolPnt;
    By[solID] = BySolPnt;

    if (m_dim == 3)
    {
      Bz[solID] = BzSolPnt;
    }

    // spherical components
    const CFreal x = coords[0];
    const CFreal y = coords[1];
    const CFreal z = (m_dim == 3) ? coords[2] : 0.0;
    const CFreal r = sqrt(x*x + y*y + z*z);
    const CFreal rXY = sqrt(x*x + y*y);

    Br[solID] = x/r*BxSolPnt + y/r*BySolPnt + z/r*BzSolPnt;
    Btheta[solID] = -y*BxSolPnt + x*BySolPnt;
    Bphi[solID] = z*x/rXY*BxSolPnt + z*y/rXY*BySolPnt - rXY*BzSolPnt;
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  DiffRHSJacobFluxReconstruction::setup();
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get the number of cells in the mesh
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  // get datahandle
  DataHandle< CFreal > Bx = socket_Bx.getDataHandle();
  DataHandle< CFreal > By = socket_By.getDataHandle();
  DataHandle< CFreal > Bz = socket_Bz.getDataHandle();
  DataHandle< CFreal > Br = socket_Br.getDataHandle();
  DataHandle< CFreal > Btheta = socket_Btheta.getDataHandle();
  DataHandle< CFreal > Bphi = socket_Bphi.getDataHandle();
  
  const CFuint nbStates = nbrCells*m_nbrSolPnts;

  // resize socket
  Bx.resize(nbStates);
  By.resize(nbStates);
  Bz.resize(nbStates);
  Br.resize(nbStates);
  Btheta.resize(nbStates);
  Bphi.resize(nbStates);
  
  // get the diffusive varset
  m_diffVarSetPoisson = m_diffusiveVarSet.d_castTo< Physics::Poisson::PoissonDiffVarSet >();
  cf_assert(m_diffVarSetPoisson.isNotNull());
  
  m_convVarSetPoisson = getMethodData().getUpdateVar().d_castTo<Physics::Poisson::PoissonConvVarSet>();
  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

