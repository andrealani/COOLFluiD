#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/ConvVolTermRHSJacobSpectralFD.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

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

MethodCommandProvider<ConvVolTermRHSJacobSpectralFD, SpectralFDMethodData, SpectralFDModule> ConvVolTermRHSJacobSpectralFDProvider("ConvVolTermRHSJacob");

//////////////////////////////////////////////////////////////////////////////

ConvVolTermRHSJacobSpectralFD::ConvVolTermRHSJacobSpectralFD(const std::string& name) :
  ConvVolTermRHSSpectralFD(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(),
  m_pertResUpdates(),
  m_derivResUpdates()

{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ConvVolTermRHSJacobSpectralFD::~ConvVolTermRHSJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  ConvVolTermRHSSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // set the volume term data for this element type
    setVolumeTermData();

    // resize the variable m_resUpdates
    resizeResAndGradUpdates();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if cell is parallel updatable or the gradients have to be computed,
      // set cell data and reconstruct states
      if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
      {
        // set the current cell and compute the cell data in the volume term computer
        m_volTermComputer->setCurrentCell(m_cell);
        m_volTermComputer->computeCellData();

        // reconstruct the solution in the flux points
        m_volTermComputer->reconstructStates(*m_cellStates);
      }

      // if cell is parallel updatable, compute the volume term
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // compute the volume term
        m_volTermComputer->computeCellConvVolumeTerm(m_resUpdates);

        // update rhs
        updateRHS();

        // compute the convective volume term contribution to the jacobian
        computeJacobConvVolTerm();
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradients();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSJacobSpectralFD::resizeResAndGradUpdates()
{
  // call function in parent class
  ConvVolTermRHSSpectralFD::resizeResAndGradUpdates();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // number of entries in residual updates
  const CFuint nbrResUpdates = m_nbrEqs*( sdLocalData[m_iElemType]->getNbrOfSolPnts() );

  // resize m_pertResUpdates and m_derivResUpdate
  m_pertResUpdates .resize(nbrResUpdates);
  m_derivResUpdates.resize(nbrResUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSJacobSpectralFD::setVolumeTermData()
{
  ConvVolTermRHSSpectralFD::setVolumeTermData();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(nbrSolPnts,nbrSolPnts,m_nbrEqs));
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSJacobSpectralFD::computeJacobConvVolTerm()
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

  // loop over the states/solpnts in this cell to perturb the states
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
      m_volTermComputer->backupAndReconstructPhysVar(iEqPert,*m_cellStates);

      // compute the perturbed volume term
      m_volTermComputer->computeCellConvVolumeTerm(m_pertResUpdates);

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

      // restore physical variable in the flux points
      m_volTermComputer->restorePhysVar(iEqPert);
    }
  }
//   acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  ConvVolTermRHSSpectralFD::setup();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  ConvVolTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
