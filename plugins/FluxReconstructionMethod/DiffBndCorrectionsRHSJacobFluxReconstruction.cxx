#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstruction, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionModule >
DiffBndCorrectionsRHSJacobFluxReconstructionProvider("DiffBndCorrectionsRHSJacob");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstruction::DiffBndCorrectionsRHSJacobFluxReconstruction(const std::string& name) :
  DiffBndCorrectionsRHSFluxReconstruction(name),
  m_cellBuilder(CFNULL),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_derivResUpdates(),
  m_pertCorrections(),
  m_resUpdates(),
  m_pertCellStatesFlxPnt(),
  m_pertSol(),
  m_pertVar(),
  m_solFlxDep(CFNULL),
  m_nbrFlxDep(),
  m_cellStatesFlxPntBackup(),
  m_influencedFlxPnts(),
  m_NbInfluencedFlxPnts(),
  elemShape(),
  m_cellStatesFlxPnt2()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstruction::~DiffBndCorrectionsRHSJacobFluxReconstruction()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  DiffBndCorrectionsRHSFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCFluxReconstruction TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  
  CFLog(VERBOSE,"DiffBndCorrectionsRHSJacobFluxReconstruction::executeOnTRS: " << faceTrs->getName() << "\n");

  // get bndFacesStartIdxs from FluxReconstructionMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;
  
  // get the geodata of the cell builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCB = m_cellBuilder->getDataGE();
  geoDataCB.trs = cellTrs;
  
  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    nbOrients = bndFacesStartIdxs[iTR].size()-1;
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
      
      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts=(*m_faceFlxPntConn)[m_orient].size();

      // select the correct flx pnts on the face out of all cell flx pnts for the current orient
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntsLocalCoords[iFlx] = (*m_allCellFlxPnts)[(*m_faceFlxPntConn)[m_orient][iFlx]];
      }
      
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();
	
	// get the states in the neighbouring cell
        m_cellStates = m_face->getNeighborGeo(0)->getStates();
	
	// compute volume
        m_cellVolume = m_face->getNeighborGeo(0)->computeVolume();
	
	cf_assert(m_cellVolume > 0.0);

        // if cell is parallel updatable, compute the correction flux
        if ((*m_cellStates)[0]->isParUpdatable())
        {
	  // build the neighbouring cell
          const CFuint cellID = m_face->getNeighborGeo(0)->getID();
          geoDataCB.idx = cellID;
          m_intCell = m_cellBuilder->buildGE();
          
          // metrics of the interior cell, needed by the compact boundary face gradient
          setIntCellMetrics();
	  
	  // set the bnd face data
	  setBndFaceData(m_face->getID());//faceID

	  // compute the states and ghost states in the flx pnts
	  computeFlxPntStates();

          // the compact BR2 face gradient of this boundary face
          computeCompactBR2BndFaceGradient();

	  // compute FI
          computeInterfaceFlxCorrection();

          // compute the wave speed updates
          computeWaveSpeedUpdates(m_waveSpeedUpd);

          // update the wave speeds
          updateWaveSpeed();

	  // compute the correction -(FI)divh of the bnd face for each sol pnt
          computeCorrection(m_corrections);

	  // update the rhs
          updateRHS();
          
          const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    
          const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
    
          const CFuint interval = iter - iterFreeze;
      
          // no perturbation loop when no Jacobian is wanted (JFNK matvecs, residual only evaluations)
          if (getMethodData().doComputeJacobian() &&
              (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0))
          {

	    // compute the contribution to the jacobian
            computeJacobDiffBndContribution();
          
          }

	  // release the cell
          m_cellBuilder->releaseGE();
        } 
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeJacobDiffBndContribution()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
    
    // put the perturbed and unperturbed corrections in the correct format
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[m_nbrEqs*iSol+iVar] = m_corrections[iSol][iVar];
    }
  }
  
  /// store values that will be overwritten
  storeBackups();

  // loop over the states in the internal cell to perturb the states
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[m_pertSol];

    // Loop over flux points to determine which flx pnts are influenced by the pert
    m_influencedFlxPnts.resize(0);
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // get current flx pnt idx
      const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
      for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
      {
        if (currFlxIdx == (*m_solFlxDep)[m_pertSol][jFlxPnt])
        {
          m_influencedFlxPnts.push_back(iFlxPnt);
        }
      }
    }
    m_NbInfluencedFlxPnts= m_influencedFlxPnts.size();

    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      // perturb physical variable in state
      m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

      // compute the perturbed states, ghost states and compact face gradient
      extrapolatePerturbedState();

      // compute the perturbed interface flx correction
      computeInterfaceFlxCorrection();
      computeCorrection(m_pertCorrections);

      // put the perturbed and unperturbed corrections in the correct format
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	{
          m_pertResUpdates[m_nbrEqs*iState+iVar] = m_pertCorrections[iState][iVar];
        }
      }

      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
      {
        acc.addValues(iSol,m_pertSol,m_pertVar,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[m_pertVar]);
      
      /// restore overwritten values
      restoreFromBackups();
    }
  }
   //acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix (or direct element blocks)
    getMethodData().assembleJacobBlock(acc, m_intCell->getID());
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::extrapolatePerturbedState()
{
  // The states in the flux points have to follow the perturbation as well as
  // the gradients: the boundary diffusive flux depends on the state at the
  // face through the transport properties and through the ghost state.
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    m_nbrSolDep = ((*m_flxSolDep)[currFlxIdx]).size();
    *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];
      *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]));
    }
  }
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);

  // rebuild the same compact boundary gradient the unperturbed residual used,
  // from the perturbed nodal states; the ghost states above are needed first,
  // the boundary value rule reads them
  computeCompactBR2BndFaceGradient();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::storeBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_cellStatesFlxPntBackup[iFlxPnt] = *(m_cellStatesFlxPnt[iFlxPnt]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::restoreFromBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  { 
    (*(m_cellStatesFlxPnt[m_influencedFlxPnts[iFlxPnt]]))[m_pertVar] = m_cellStatesFlxPntBackup[m_influencedFlxPnts[iFlxPnt]][m_pertVar];  
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTerm(iEq,iFlx) = (*(m_cellStatesFlxPnt[iFlx]->getData()))[iEq];
      ghostGradTerm(iEq,iFlx) = (*(m_flxPntGhostSol[iFlx]->getData()))[iEq];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTerm(iEq,iFlx) = (*(m_cellStatesFlxPnt2[iFlx]->getData()))[iEq];
      ghostGradTerm(iEq,iFlx) = (*(m_flxPntGhostSol[iFlx]->getData()))[iEq];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeCellGradTerm(RealMatrix& gradTerm)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iSol = 0; iSol < m_nbrFaceFlxPnts; ++iSol)
    {
      gradTerm(iEq,iSol) = (*((*m_cellStates)[iSol]->getData()))[iEq];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTermL(iEq,iFlx) = (*(m_pertCellStatesFlxPnt[LEFT][iFlx]->getData()))[iEq];
      gradTermR(iEq,iFlx) = (*(m_pertCellStatesFlxPnt[RIGHT][iFlx]->getData()))[iEq];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeFlxPntGradTerm2(RealMatrix& gradTerm)
{
  // the base class reconstructs the update variables themselves, so the
  // transform is the identity here
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTerm(iEq,iFlx) = (*(m_cellStatesFlxPnt2[iFlx]->getData()))[iEq];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  DiffBndCorrectionsRHSFluxReconstruction::setup();
  
  // get CellToFaceGeBuilder
  m_cellBuilder      = getMethodData().getCellBuilder();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);

  //get element shape
  elemShape = frLocalData[0]->getShape();
  
  m_solFlxDep = frLocalData[0]->getSolPntFlxDependency();

  m_nbrFlxDep = ((*m_solFlxDep)[0]).size();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));

  // resize variables
  const CFuint nbrCellResiduals = m_nbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(nbrCellResiduals);
  m_derivResUpdates.resize(nbrCellResiduals);
  m_resUpdates.resize(nbrCellResiduals);
  m_pertCorrections.resize(m_nbrSolPnts);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pertCorrections[iSol].resize(m_nbrEqs);
  }
  
  // create internal and ghost states
  m_pertCellStatesFlxPnt.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_pertCellStatesFlxPnt[LEFT].push_back(new State());
    m_pertCellStatesFlxPnt[RIGHT].push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_pertCellStatesFlxPnt[LEFT][iFlx]->setLocalID(iFlx);
    m_pertCellStatesFlxPnt[RIGHT][iFlx]->setLocalID(iFlx);
  }
  
  m_cellStatesFlxPntBackup.resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPntBackup[iFlx].resize(m_nbrEqs); 
  }
  
  // create internal states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt2.push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt2[iFlx]->setLocalID(iFlx);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt2[iFlx]);
    deletePtr(m_pertCellStatesFlxPnt[LEFT][iFlx]);
    deletePtr(m_pertCellStatesFlxPnt[RIGHT][iFlx]);
  }
  m_pertCellStatesFlxPnt.clear();
  m_cellStatesFlxPnt2.clear();

  // unsetup parent class
  DiffBndCorrectionsRHSFluxReconstruction::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
