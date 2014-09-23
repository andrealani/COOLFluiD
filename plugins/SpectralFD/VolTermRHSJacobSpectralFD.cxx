#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/VolTermRHSJacobSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<VolTermRHSJacobSpectralFD, SpectralFDMethodData, SpectralFDModule> VolTermRHSJacobSpectralFDProvider("VolTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

VolTermRHSJacobSpectralFD::VolTermRHSJacobSpectralFD(const std::string& name) :
  VolTermRHSSpectralFD(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(),
  m_pertResUpdates(),
  m_derivResUpdates(),
  m_gradsMBndFaceTerm(),
  m_physVarGradUpdates()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

VolTermRHSJacobSpectralFD::~VolTermRHSJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  VolTermRHSSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // SET THE VOLUME AND FACE TERM COMPUTERS DATA FOR THIS ELEMENT TYPE
    setVolumeAndFaceTermComputersData();

    // RESIZE THE VARIABLES M_RESUPDATES AND M_GRADUPDATES
    resizeResAndGradUpdates();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if cell is parallel updatable, compute the volume term
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // COMPUTE DATA IN VOLUME AND FACE TERM COMPUTERS
        computeCellVolumeAndFaceTermData();

        // COMPUTE CELL GRADIENTS
        computeAndReconstructGradients();

        // COMPUTE CONVECTIVE VOLUME TERMS
        m_volTermComputer->computeCellConvVolumeTerm(m_resUpdates);

        // COMPUTE DIFFUSIVE VOLUME TERMS
        m_volTermComputer->computeCellDiffVolumeTerm(m_diffResUpdates);

        // ADD TOTAL UPDATE TO RESIDUAL
        m_resUpdates += m_diffResUpdates;
        addUpdatesToResidual();

        // COMPUTE JACOBIAN CONTRIBUTION OF VOLUME TERMS
        computeCellGradsMinusBndFaceTerms();
        computeJacobVolTerm();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  CFTRACEEND;
 }

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::resizeResAndGradUpdates()
{
  // call parent class function
  VolTermRHSSpectralFD::resizeResAndGradUpdates();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // number of entries in residual updates
  const CFuint nbrResUpdates = m_nbrEqs*nbrSolPnts;

  // resize m_pertResUpdates and m_derivResUpdate
  m_pertResUpdates .resize(nbrResUpdates);
  m_derivResUpdates.resize(nbrResUpdates);

  // resize m_physVarGradUpdates
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_physVarGradUpdates[iSide].resize(nbrSolPnts);
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      m_physVarGradUpdates[iSide][iSol].resize(m_dim);
    }
  }

  // resize m_gradsMBndFaceTerm
  m_gradsMBndFaceTerm.resize(nbrSolPnts);
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_gradsMBndFaceTerm[iSol].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsMBndFaceTerm[iSol][iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::setVolumeAndFaceTermComputersData()
{
  // call parent class function
  VolTermRHSSpectralFD::setVolumeAndFaceTermComputersData();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(nbrSolPnts,nbrSolPnts,m_nbrEqs));
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::computeCellGradsMinusBndFaceTerms()
{
  // copy gradients
  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsMBndFaceTerm[iSol][iEq] = (*m_cellGrads[iSol])[iEq];
    }
  }

  // subtract boundary face terms
  // compute face contributions to the gradients
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // subtract previous gradient contribution from other variable gradients
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_gradsMBndFaceTerm[iSol][iEq] -= m_bndFaceGradUpdates[iFace][iSol][iEq]/m_solJacobDet[iSol];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::computeAndReconstructPertPhysVarGrad(const CFuint iVar)
{
  // copy unperturbed gradients minus boundary face terms to perturbed gradients variable and
  // compute volume contribution to physical variable gradient
  m_volTermComputer->computePhysVarGradientVolumeTerm(iVar,m_physVarGradUpdates[0]);
  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // copy unperturbed gradients minus boundary face term contributions
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[iSol])[iEq] = m_gradsMBndFaceTerm[iSol][iEq];
    }

    // volume term of current physical variable gradient
    (*m_cellGrads[iSol])[iVar] = m_physVarGradUpdates[0][iSol];
  }

  // compute face contributions to the gradients
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // compute the boundary face contribution to the gradients
      m_bndFaceTermComputers[iFace]->computeGradientFaceTerm(m_gradUpdates[0]);

      // add the contribution to the current physical variable gradient
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        (*m_cellGrads[iSol])[iVar] += m_gradUpdates[0][iSol][iVar];
      }

      // add new gradient contribution to other variable gradients
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        if (iEq != iVar)
        {
          for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
          {
            (*m_cellGrads[iSol])[iEq] += m_gradUpdates[0][iSol][iEq]/m_solJacobDet[iSol];
          }
        }
      }
    }
    else
    {
      // compute the face contribution to the gradients
      m_faceTermComputers[iFace]->computePhysVarGradFaceTerm(iVar,m_physVarGradUpdates);

      // add the contribution to the current physical variable gradient
      const CFuint side = (*m_currCellSide)[iFace];
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        (*m_cellGrads[iSol])[iVar] += m_physVarGradUpdates[side][iSol];
      }
    }
  }

  // divide current physical variable gradient by solution point Jacobian determinant
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solJacobDet[iSol];
    (*m_cellGrads[iSol])[iVar] *= invJacobDet;
  }

  // reconstruct the gradients in the flux points
  // here, all the gradients are still always reconstructed,
  // since all of them could change through the boundary face term
  m_volTermComputer->reconstructGradients(m_cellGrads);
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::backupAndReconstructVolumeAndFacePhysVars(const CFuint iVar)
{
  // volume term computer
  m_volTermComputer->backupAndReconstructPhysVar(iVar,*m_cellStates);

  // face term computers
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // reconstruct the states in the face flux points
      m_bndFaceTermComputers[iFace]->backupAndReconstructPhysVar(iVar,*m_cellStates);
    }
    else
    {
      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->backupAndReconstructPhysVar((*m_currCellSide)[iFace],iVar,
                                                              *m_cellStates);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::restoreVolumeAndFacePhysVars(const CFuint iVar)
{
  // volume term computer
  m_volTermComputer->restorePhysVar(iVar);

  // face term computers
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if ((*m_isFaceOnBoundary)[iFace])
    {
      // reconstruct the states in the face flux points
      m_bndFaceTermComputers[iFace]->restorePhysVar(iVar);
    }
    else
    {
      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->restorePhysVar((*m_currCellSide)[iFace],iVar);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::computeJacobVolTerm()
{
  // get number of solution points in this cell
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

  // loop over the states in this cell to perturb the states
  for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[iSolPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
    {
      // perturb physical variable in state
      m_numJacob->perturb(iEqPert,pertState[iEqPert]);

      // backup and reconstruct physical variable in the volume and face flux points
      backupAndReconstructVolumeAndFacePhysVars(iEqPert);

      // compute and reconstruct the perturbed cell gradients after the perturbation
      computeAndReconstructPertPhysVarGrad(iEqPert);

      // compute the perturbed volume term
      m_volTermComputer->computeCellConvVolumeTerm(m_pertResUpdates);
      m_volTermComputer->computeCellDiffVolumeTerm(m_diffResUpdates);
      m_pertResUpdates += m_diffResUpdates;

      // compute the finite difference derivative of the volume term
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

      // restore physical variable in the volume and face flux points
      restoreVolumeAndFacePhysVars(iEqPert);
    }
  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  VolTermRHSSpectralFD::setup();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // resize m_physVarGradUpdates
  m_physVarGradUpdates.resize(2);
}

//////////////////////////////////////////////////////////////////////////////

void VolTermRHSJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of the parent class
  VolTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
