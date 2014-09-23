#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/ConvVolTermDiagBlockJacobSpectralFD.hh"
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

MethodCommandProvider<ConvVolTermDiagBlockJacobSpectralFD, SpectralFDMethodData, SpectralFDModule> ConvVolTermDiagBlockJacobSpectralFDProvider("ConvVolTermDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

ConvVolTermDiagBlockJacobSpectralFD::ConvVolTermDiagBlockJacobSpectralFD(const std::string& name) :
  ConvVolTermRHSJacobSpectralFD(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  m_currDiagMatrix(CFNULL)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ConvVolTermDiagBlockJacobSpectralFD::~ConvVolTermDiagBlockJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermDiagBlockJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermDiagBlockJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  ConvVolTermRHSJacobSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermDiagBlockJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

  // get the diagonal block Jacobian matrix datahandle
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

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
        // set current diagonal block Jacobian matrix
        m_currDiagMatrix = &diagBlockJacobMatr[elemIdx];

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

        // compute the convective volume term contribution to the diagonal block jacobian
        computeDiagBlockJacobConvVolTerm();
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

void ConvVolTermDiagBlockJacobSpectralFD::setVolumeTermData()
{
  // call the setVolumeTermData() of the parent class of the parent class,
  // such that the accumulator for the lss matrix is not set.
  ConvVolTermRHSSpectralFD::setVolumeTermData();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermDiagBlockJacobSpectralFD::computeDiagBlockJacobConvVolTerm()
{
  // get number of solution points in this cell
  const CFuint nbrSolPnts = m_cellStates->size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // loop over the states/solpnts in this cell to perturb the states
  CFuint pertResUpdIdx = 0;
  for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[iSolPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert, ++pertResUpdIdx)
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

      // add the derivative of the residual updates to diagonal block matrix
      addToDiagBlockJacobMatrix(pertResUpdIdx);

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the flux points
      m_volTermComputer->restorePhysVar(iEqPert);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermDiagBlockJacobSpectralFD::addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx)
{
  const CFuint nbrRes = m_currDiagMatrix->nbRows();
  cf_assert(nbrRes == m_derivResUpdates.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    (*m_currDiagMatrix)(iRes,pertResUpdIdx) += m_derivResUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermDiagBlockJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // call setup of parent class of parent class, in order not to get the lss
  ConvVolTermRHSSpectralFD::setup();

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermDiagBlockJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of parent class of parent class, in order not to get the lss
  ConvVolTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    ConvVolTermDiagBlockJacobSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ConvVolTermRHSJacobSpectralFD::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr     );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
