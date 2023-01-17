#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstruction.hh"
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

MethodCommandProvider< ConvBndCorrectionsRHSJacobFluxReconstruction, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvBndCorrectionsRHSJacobFluxReconstructionProvider("ConvBndCorrectionsRHSJacob");

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSJacobFluxReconstruction::ConvBndCorrectionsRHSJacobFluxReconstruction(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstruction(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_resUpdates(),
  m_derivResUpdates(),
  m_pertCorrections(),
  m_pertSol(),
  m_pertVar(),
  m_solFlxDep(CFNULL),
  m_nbrFlxDep(),
  m_cellStatesFlxPntBackup(),
  m_influencedFlxPnt(),
  m_flxPntRiemannFluxBackup(),
  elemShape()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSJacobFluxReconstruction::~ConvBndCorrectionsRHSJacobFluxReconstruction()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  ConvBndCorrectionsRHSFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current bnd face TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  
  CFLog(VERBOSE,"ConvBndCorrectionRHSJacobFluxReconstruction::executeOnTRS: " << faceTrs->getName() << "\n");

  // get bndFacesStartIdxs from FluxReconstructionMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;
  
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();
  m_bcStateComputer->preProcess();
  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
      
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

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);
	
	// get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();
	
        CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");
	CFLog(VERBOSE,"coord state 0: " << (((*m_cellStates)[0])->getCoordinates()) << "\n");

        // if cell is parallel updatable or the gradients have to be computed, compute the necessary data
        if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
        {  
	  // set the bnd face data
	  setBndFaceData(m_face->getID());//faceID
	  
	  // compute the perturbed states and ghost states in the flx pnts
          computeFlxPntStates();
	}
	
	// if the cell is parallel updatable, compute the flx correction
	if ((*m_cellStates)[0]->isParUpdatable())
	{
	  // compute FI-FD
          computeInterfaceFlxCorrection();
	  
          // compute the wave speed updates
          computeWaveSpeedUpdates(m_waveSpeedUpd);
      
          // update the wave speeds
          updateWaveSpeed();
       
	  // compute the correction -(FI-FD)divh of the bnd face for each sol pnt
          computeCorrection(m_corrections);
	  
	  // update the rhs
          updateRHS();
	}
	  
	// if there is a diffusive term, compute the gradients
        if (hasDiffTerm)
        {
          computeGradientBndFaceCorrections();
        }
        
        const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    
        const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
    
        const CFuint interval = iter - iterFreeze;
      
        if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
        {
	
	  // if the cell is parallel updatable, compute the contribution to the numerical jacobian
	  if ((*m_cellStates)[0]->isParUpdatable())
	  {
	    // compute the convective boundary flux correction contribution to the jacobian
	    computeJacobConvBndCorrection();
          }
        }
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::computeJacobConvBndCorrection()
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
  
  // store backups of values that will be overridden
  storeBackups();

  // loop over the states in the internal cell to perturb the states
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[m_pertSol];
    
    // Loop over flux points to determine which flx pnts are influenced by the pert
    if (elemShape == CFGeoShape::TRIAG || elemShape == CFGeoShape::TETRA)
    {
      m_influencedFlxPnt = 0;
      m_NbInfluencedFlxPnts = m_nbrFaceFlxPnts;
    }
    else
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get current flx pnt idx
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
      
        for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
        {
          if (currFlxIdx == (*m_solFlxDep)[m_pertSol][jFlxPnt])
          {
            m_influencedFlxPnt = iFlxPnt;
            m_NbInfluencedFlxPnts= iFlxPnt+1;
          }
        }
      }
    }

    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      // perturb physical variable in state
      m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
      
      // compute the perturbed states and ghost states in the flx pnts
      extrapolatePerturbedState();

      // compute the perturbed interface flx correction
      computePertInterfaceFlxCorrection();
      computePertCorrection(m_pertCorrections);
      
      // put the perturbed and unperturbed corrections in the correct format
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	{
          m_pertResUpdates[m_nbrEqs*iState+iVar] = m_pertCorrections[iState][iVar];
        }
      }

      // add contributions to the Jacobian
      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;
      
      for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
      {   
        const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
          
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];
          acc.addValues(solIdx,m_pertSol,m_pertVar,&m_derivResUpdates[m_nbrEqs*solIdx]);
        }  
      } 

      // restore physical variable in state
      m_numJacob->restore(pertState[m_pertVar]);
      
      // restore the overridden values
      restoreFromBackups();
    }
  }
  
//   if (m_intCell->getID() == 56)
//   {
//   CFLog(VERBOSE,"accBndFace:\n");
//    acc.printToScreen();
//   }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::extrapolatePerturbedState()
{       
  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {   
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    // reset the extrapolated states
    (*(m_cellStatesFlxPnt[iFlxPnt]))[m_pertVar] = 0.0;
    
    // extrapolate the states to current flx pnt
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

    (*(m_cellStatesFlxPnt[iFlxPnt]))[m_pertVar] += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]))[m_pertVar];
    }
  }  
  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::computePertInterfaceFlxCorrection()
{ 
  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt){   
    // compute the riemann flux in the flx pnts
    m_flxPntRiemannFlux[iFlxPnt] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[iFlxPnt]),
                          *(m_flxPntGhostSol[iFlxPnt]),
                          m_unitNormalFlxPnts[iFlxPnt]);
      
    // store the local Riemann flux, scaled with geometrical Jacobian
    m_flxPntRiemannFlux[iFlxPnt] *= m_faceJacobVecSizeFlxPnts[iFlxPnt];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::computePertCorrection(vector< RealVector >& corrections)
{ 
  cf_assert(corrections.size() == m_nbrSolPnts);  
  // reset corrections
  const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][m_influencedFlxPnt];
        
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
  {
    const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];
    corrections[solIdx] = 0.0;
  }  

  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  { 
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    // the current correction factor previously computed
    const RealVector& currentCorrFactor = m_flxPntRiemannFlux[iFlxPnt];

    cf_assert(currentCorrFactor.size() == m_nbrEqs);

    // compute the term due to each flx pnt
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

      // divergence of the correctionfct
      const CFreal divh = m_corrFctDiv[solIdx][flxIdx];

      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {  
        corrections[solIdx][iVar] -= currentCorrFactor[iVar] * divh; 
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::storeBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_cellStatesFlxPntBackup[iFlxPnt] = *(m_cellStatesFlxPnt[iFlxPnt]);
    m_flxPntRiemannFluxBackup[iFlxPnt] = m_flxPntRiemannFlux[iFlxPnt];
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pertCorrections[iSol] = m_corrections[iSol];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::restoreFromBackups()
{
  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  { 
  (*(m_cellStatesFlxPnt[iFlxPnt]))[m_pertVar] = m_cellStatesFlxPntBackup[iFlxPnt][m_pertVar];
  m_flxPntRiemannFlux[iFlxPnt] = m_flxPntRiemannFluxBackup[iFlxPnt];
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
  {
    m_pertCorrections[iSol] = m_corrections[iSol];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  ConvBndCorrectionsRHSFluxReconstruction::setup();
  
  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  //get element shape
  elemShape = frLocalData[0]->getShape();

  m_solFlxDep = frLocalData[0]->getSolPntFlxDependency();

  m_nbrFlxDep = ((*m_solFlxDep)[0]).size();
  
  // resize variables
  const CFuint nbrRes = m_nbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(nbrRes);
  m_derivResUpdates.resize(nbrRes);
  m_resUpdates .resize(nbrRes);
  m_pertCorrections.resize(m_nbrSolPnts);
  m_cellStatesFlxPntBackup.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxBackup.resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPntBackup[iFlx].resize(m_nbrEqs); 
    m_flxPntRiemannFluxBackup[iFlx].resize(m_nbrEqs);
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pertCorrections[iSol].resize(m_nbrEqs);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  ConvBndCorrectionsRHSFluxReconstruction::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
