// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionPoisson/ConvDiffJacobFluxReconstructionPoisson.hh"
#include "FluxReconstructionPoisson/FluxReconstructionPoisson.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "Poisson/PoissonDiffVarSet.hh"
#include "Poisson/PoissonConvVarSet.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::Poisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvDiffJacobFluxReconstructionPoisson,
		       FluxReconstructionSolverData,
		       FluxReconstructionPoissonModule >
convDiffRHSJacobPoissonFluxReconstructionProvider("ConvDiffRHSJacobPoisson");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffJacobFluxReconstructionPoisson::ConvDiffJacobFluxReconstructionPoisson(const std::string& name) :
  ConvDiffJacobFluxReconstruction(name),
  socket_Bx("Bx"),
  socket_By("By"),
  socket_Bz("Bz"),
  socket_Br("Br"),
  socket_Btheta("Btheta"),
  socket_Bphi("Bphi"),
  m_pData()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionPoisson::configure ( Config::ConfigArgs& args )
{
  ConvDiffJacobFluxReconstruction::configure(args);
} 

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  ConvDiffJacobFluxReconstructionPoisson::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = ConvDiffJacobFluxReconstruction::providesSockets();

  result.push_back(&socket_Bx);
  result.push_back(&socket_By);
  result.push_back(&socket_Bz);
  result.push_back(&socket_Br);
  result.push_back(&socket_Btheta);
  result.push_back(&socket_Bphi);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionPoisson::computeInterfaceFlxCorrection()
{
  // common diffusive flux with the compact face gradients, no convective Riemann flux is subtracted
  DiffRHSFluxReconstruction::computeInterfaceFlxCorrection();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntRiemannFluxDiff[iFlx] = m_flxPntRiemannFlux[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionPoisson::prepareFluxComputation()
{
  //const bool isPerturb = this->getMethodData().isPerturb();
  //const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  //m_diffVarSetPoisson->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionPoisson::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // here convective and artificial parts are added!
  //ConvDiffJacobFluxReconstruction::computeWaveSpeedUpdates(waveSpeedUpd);
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
  }
          
  // now add diffusive part
  CFreal visc = 33.0;

  //const CFreal dynVisc = m_navierStokesVarSetPoisson->getCurrDynViscosity();
  
  //const CFreal factorPr = min(m_navierStokesVarSetPoisson->getModel().getPrandtl(),1.0);
  //cf_assert(factorPr>0.0);
    
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    //for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      //const CFreal rho = m_navierStokesVarSetPoisson->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      //visc = dynVisc/rho/factorPr;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionPoisson::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // side that is not perturbed
    const CFuint iOtherSide = (m_pertSide == LEFT) ? RIGHT : LEFT;

    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      // state extrapolated to the flux point on the perturbed side
      State& pertState = *(m_cellStatesFlxPnt[m_pertSide][iFlx]);

      // average of the compact face gradients of the two sides
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        *(m_avgGrad[iEq]) = (*(m_cellGradFlxPnt[LEFT][iFlx][iEq]) + *(m_cellGradFlxPnt[RIGHT][iFlx][iEq]))/2.0;
      }

      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

        // average of the perturbed state and the state of the other side
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_avgSol[iEq] = (pertState[iEq] + (*(m_cellStatesFlxPnt[iOtherSide][iFlx]))[iEq])/2.0;
        }

        prepareFluxComputation();

        // perturbed common diffusive flux
        computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlx],0,m_flxPntRiemannFluxPert[iFlx]);

        // derivative of the common flux times the residual factor
        m_numJacob->computeDerivative(m_flxPntRiemannFluxPert[iFlx],m_flxPntRiemannFlux[iFlx],m_riemannFluxJacobian[m_pertSide][iFlx][m_pertVar]);
        m_riemannFluxJacobian[m_pertSide][iFlx][m_pertVar] *= resFactor;

        m_numJacob->restore(pertState[m_pertVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionPoisson::computeUnpertCellDiffResiduals(const CFuint side)
{
  // diffusive residual of the cell
  ConvDiffJacobFluxReconstruction::computeUnpertCellDiffResiduals(side);

  // magnetic field of the cell
  computeMagneticField(side);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionPoisson::computeMagneticField(const CFuint side)
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

void ConvDiffJacobFluxReconstructionPoisson::setup()
{
  CFAUTOTRACE;

  // setup parent class
  ConvDiffJacobFluxReconstruction::setup();
  
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

  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();

  m_updateToSolutionVecTrans->setup(2);

  // get Euler varset
  m_convVarSetPoisson = getMethodData().getUpdateVar().d_castTo<Physics::Poisson::PoissonConvVarSet>();

  m_convVarSetPoisson->getModel()->resizePhysicalData(m_pData);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
