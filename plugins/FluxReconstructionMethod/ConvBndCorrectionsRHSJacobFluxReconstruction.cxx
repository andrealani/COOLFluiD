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
  m_pertCorrections()
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

  // get current QuadFreeBCFluxReconstruction TRS
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
	
	CFLog(VERBOSE,"faceID: " << faceID << "\n");
	CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");
	CFLog(VERBOSE,"coord state 0: " << (((*m_cellStates)[0])->getCoordinates()) << "\n");
	
	setBndFaceData(m_face->getID());//faceID

        // if cell is parallel updatable, compute the correction flux
        if ((*m_cellStates)[0]->isParUpdatable())
        {
	  // compute the discontinuous flx in the flx pnts
	  computeDiscontinuousFlx();
	  
	  // set the face ID in the BCStateComputer
	  m_bcStateComputer->setFaceID(m_face->getID());
      
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
	  
	  // compute the convective boundary flux correction contribution to the jacobian
	  computeJacobConvBndCorrection();
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
  // get number of solution points in the cell
  const CFuint nbrSolPnts = m_cellStates->size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
  }
  
  // put the perturbed and unperturbed corrections in the correct format
  for (CFuint iState = 0; iState < nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[m_nbrEqs*iState+iVar] = m_corrections[iState][iVar];
    }
  }

  // loop over the states in the internal cell to perturb the states
  for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[iSolPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
    {
      // perturb physical variable in state
      m_numJacob->perturb(iEqPert,pertState[iEqPert]);

      // backup and reconstruct physical variable in the flux points
      // and reconstruct the ghost states)
      backupAndReconstructPhysVar(iEqPert,*m_cellStates);
      
      computeDiscontinuousFlx();

      // compute the perturbed boundary face term
      computeInterfaceFlxCorrection();
      computeCorrection(m_pertCorrections);
      
      // get nb of states
      const CFuint nbrStates = m_cellStates->size();
      
      // put the perturbed and unperturbed corrections in the correct format
      for (CFuint iState = 0; iState < nbrStates; ++iState)
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

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
      {
        acc.addValues(iSol,iSolPert,iEqPert,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the flux points
      restorePhysVar(iEqPert);
      
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

void ConvBndCorrectionsRHSJacobFluxReconstruction::backupAndReconstructPhysVar
    (const CFuint iVar, const vector< State* >& cellIntStates)
{
  // backup
  backupPhysVar(iVar);
  

}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::backupPhysVar(const CFuint iVar)
{
  //cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVar.size());
  //cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStates.size());
  //for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  //{
    //m_backupPhysVar[iFlx] = (*m_flxPntIntSol[iFlx])[iVar];
    //m_backupGhostStates[iFlx] = *m_flxPntGhostSol[iFlx];
  //}
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::restorePhysVar(const CFuint iVar)
{
  //cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVar.size());
  //cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStates.size());
  //for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  //{
    //(*m_flxPntIntSol[iFlx])[iVar] = m_backupPhysVar[iFlx];
    //*m_flxPntGhostSol[iFlx] = m_backupGhostStates[iFlx];
  //}
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

  ConvBndCorrectionsRHSFluxReconstruction::setup();
  
  // get the states reconstructor
  m_statesReconstr = getMethodData().getStatesReconstructor();
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  
  // number of sol points
  const CFuint nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
  
  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));
  
  // resize variables
  const CFuint nbrRes = m_nbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(nbrRes);
  m_derivResUpdates.resize(nbrRes);
  m_resUpdates .resize(nbrRes);
  m_pertCorrections.resize(nbrSolPnts);
  
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_pertCorrections[iSol].resize(m_nbrEqs);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;

  ConvBndCorrectionsRHSFluxReconstruction::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
