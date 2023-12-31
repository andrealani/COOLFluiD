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
  m_dampCoeff(),
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

void DiffRHSJacobFluxReconstructionPoisson::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "DiffRHSJacobFluxReconstructionPoisson::execute()\n");
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;
  
  // reset cell flags
  for (CFuint iCell = 0; iCell < m_cellFlags.size(); ++iCell)
  {
    m_cellFlags[iCell] = false;
  }
  
  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;
  
  // get the geodata of the cell builders and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCBL = m_cellBuilders[LEFT]->getDataGE();
  geoDataCBL.trs = cells;
  CellToFaceGEBuilder::GeoData& geoDataCBR = m_cellBuilders[RIGHT]->getDataGE();
  geoDataCBR.trs = cells;
  
  //// Loop over faces to calculate fluxes and interface fluxes in the flux points
  
  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    CFLog(VERBOSE, "Orient = " << m_orient << "\n");
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
    m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();
      
      // compute volume
      m_cellVolume[LEFT] = m_cells[LEFT]->computeVolume();
      m_cellVolume[RIGHT] = m_cells[RIGHT]->computeVolume();
      
      cf_assert(m_cellVolume[LEFT] > 0.0);
      cf_assert(m_cellVolume[RIGHT] > 0.0);
      
      // if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
	// build the neighbouring cells
        const CFuint cellIDL = m_face->getNeighborGeo(LEFT)->getID();
        geoDataCBL.idx = cellIDL;
        m_cells[LEFT] = m_cellBuilders[LEFT ]->buildGE();
        const CFuint cellIDR = m_face->getNeighborGeo(RIGHT)->getID();
        geoDataCBR.idx = cellIDR;
        m_cells[RIGHT] = m_cellBuilders[RIGHT]->buildGE();

	// set the face data
	setFaceData(m_face->getID());//faceID

	// compute the left and right states and gradients in the flx pnts
	computeFlxPntStatesAndGrads();

	// compute FI
	computeInterfaceFlxCorrection();

	// compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);

        // update the wave speed
        updateWaveSpeed();

	// compute the correction for the left neighbour
	computeCorrection(LEFT, m_divContFlxL);
	m_divContFlx = m_divContFlxL;

	// update RHS
	updateRHS();
	
	// compute the correction for the right neighbour
	computeCorrection(RIGHT, m_divContFlxR);
	m_divContFlx = m_divContFlxR;
	
	// update RHS
	updateRHS();
        
        // compute needed cell contributions
        if (!m_cellFlags[cellIDL])
        {
          computeUnpertCellDiffResiduals(LEFT);
          m_unpertAllCellDiffRes[cellIDL] = m_unpertCellDiffRes[LEFT];

	  // update RHS
	  if ((*m_states[LEFT ])[0]->isParUpdatable()) updateRHSUnpertCell(LEFT);
        }
        if (!m_cellFlags[cellIDR])
        {
          computeUnpertCellDiffResiduals(RIGHT);
          m_unpertAllCellDiffRes[cellIDR] = m_unpertCellDiffRes[RIGHT];

	  // update RHS
	  if ((*m_states[RIGHT])[0]->isParUpdatable()) updateRHSUnpertCell(RIGHT);
        }

	// get all the faces neighbouring the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

        // set the local indexes of the other faces than the current faces
        setOtherFacesLocalIdxs();

	// make a back up of the grads and put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_cellGradsBackUp[LEFT][iState][iVar] = (*m_cellGrads[LEFT][iState])[iVar];
            m_cellGradsBackUp[RIGHT][iState][iVar] = (*m_cellGrads[RIGHT][iState])[iVar];
            
            m_resUpdates[LEFT][m_nbrEqs*iState+iVar] = m_divContFlxL[iState][iVar];
            m_resUpdates[RIGHT][m_nbrEqs*iState+iVar] = m_divContFlxR[iState][iVar];
          }
        }

	for (CFuint iSide = 0; iSide < 2; ++iSide)
        {
          // compute solution points Jacobian determinants
          m_solJacobDet[iSide] = m_cells[iSide]->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
	}

        // compute the diffusive face term contribution to the jacobian
        if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          computeBothJacobsDiffFaceTerm();
        }
        else if ((*m_states[LEFT])[0]->isParUpdatable())
        {
          computeOneJacobDiffFaceTerm(LEFT );
        }
        else if ((*m_states[RIGHT])[0]->isParUpdatable())
        {
          computeOneJacobDiffFaceTerm(RIGHT);
        }

        // release the cells
        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();
        
        m_cellFlags[cellIDL] = true;
        m_cellFlags[cellIDR] = true;
      }
      
      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
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

void DiffRHSJacobFluxReconstructionPoisson::computeInterfaceFlxCorrection()
{
    
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffVarSetPoisson->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffVarSetPoisson->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts);
    
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
       
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }

    // damping factor
    const CFreal dampFactor = m_dampCoeff*m_faceInvCharLengths[iFlxPnt];

    // compute averaged (damped) gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
      *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
    }
    
    prepareFluxComputation();
     
    // compute the Riemann flux
//     m_flxPntRiemannFlux[iFlxPnt] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0);
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  vector< RealVector* > tempStates;
  vector< RealVector* > tempGhostStates;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates.push_back(m_cellStatesFlxPnt[0][iFlx]->getData());
    tempGhostStates.push_back(m_flxPntGhostSol[iFlx]->getData());
  }
  
  m_diffVarSetPoisson->setGradientVars(tempStates,gradTerm,m_nbrFaceFlxPnts);
  m_diffVarSetPoisson->setGradientVars(tempGhostStates,ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::computeCellGradTerm(RealMatrix& gradTerm)
{ 
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  
  vector< RealVector* > tempStates;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    tempStates.push_back((*m_cellStates)[iSol]->getData());
  }
  
  m_diffVarSetPoisson->setGradientVars(tempStates,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  vector< vector< RealVector* > > tempStates;
  tempStates.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates[LEFT].push_back(m_cellStatesFlxPnt[LEFT][iFlx]->getData());
    tempStates[RIGHT].push_back(m_cellStatesFlxPnt[RIGHT][iFlx]->getData());
  }
  
  m_diffVarSetPoisson->setGradientVars(tempStates[LEFT],gradTermL,m_nbrFaceFlxPnts);
  m_diffVarSetPoisson->setGradientVars(tempStates[RIGHT],gradTermR,m_nbrFaceFlxPnts);
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
  // create a list of the dimensions in which the deriv will be calculated
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_cellFluxProjVects[iDim] = m_cells[side]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
  }
  
  // get datahandle
  DataHandle< CFreal > Bx = socket_Bx.getDataHandle();
  DataHandle< CFreal > By = socket_By.getDataHandle();
  DataHandle< CFreal > Bz = socket_Bz.getDataHandle();
  DataHandle< CFreal > Br = socket_Br.getDataHandle();
  DataHandle< CFreal > Btheta = socket_Btheta.getDataHandle();
  DataHandle< CFreal > Bphi = socket_Bphi.getDataHandle();

  // set the states
  *m_cellStates = *(m_states[side]);

  // make a backup of the grads if necessary
  if (side == RIGHT)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellGradsBackUp[0][iSol] = *(m_cellGrads[0][iSol]);
      *(m_cellGrads[0][iSol]) = *(m_cellGrads[1][iSol]);
    }
  }

  // compute the volume term
  computeDivDiscontFlx(m_pertDivContFlx[0]);

  // put the unpert discontinuous diff residual in the correct format
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_unpertCellDiffRes[side][m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
    }
    
    const CFuint stateLocalID = (((*(m_states[side]))[iState]))->getLocalID();
    
    Bx[stateLocalID] = (*(m_cellGrads[0][iState]))[0][0];
    By[stateLocalID] = (*(m_cellGrads[0][iState]))[0][1];
    if (m_dim == 3) Bz[stateLocalID] = (*(m_cellGrads[0][iState]))[0][2];
    
    const CFreal x = ((((*(m_states[side]))[iState]))->getCoordinates())[0];
    const CFreal y = ((((*(m_states[side]))[iState]))->getCoordinates())[1];
    const CFreal z = (m_dim == 3) ? ((((*(m_states[side]))[iState]))->getCoordinates())[2] : 0.0;
    const CFreal r = sqrt(x*x+y*y+z*z);
    const CFreal r2D = sqrt(x*x+y*y);

    const CFreal Bxv = (*(m_cellGrads[0][iState]))[0][0];
    const CFreal Byv = (*(m_cellGrads[0][iState]))[0][1];
    const CFreal Bzv = (m_dim == 3) ? (*(m_cellGrads[0][iState]))[0][2] : 0.0; 
    
    Br[stateLocalID] = x/r*Bxv + y/r*Byv + z/r*Bzv;
    
    Btheta[stateLocalID] = -y*Bxv + x*Byv;
    
    Bphi[stateLocalID] = z*x/r2D*Bxv + z*y/r2D*Byv - r2D*Bzv;
  }

  // restore grads if necessary
  if (side == RIGHT)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      *(m_cellGrads[0][iSol]) = m_cellGradsBackUp[0][iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionPoisson::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  DiffRHSJacobFluxReconstruction::setup();
  
  // get damping coeff
  m_dampCoeff = getMethodData().getDiffDampCoefficient();
  
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

