#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/VolTermDiagBlockJacobSpectralFD.hh"

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

MethodCommandProvider<VolTermDiagBlockJacobSpectralFD, SpectralFDMethodData, SpectralFDModule> VolTermDiagBlockJacobSpectralFDProvider("VolTermDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

VolTermDiagBlockJacobSpectralFD::VolTermDiagBlockJacobSpectralFD(const std::string& name) :
  VolTermRHSJacobSpectralFD(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  m_currDiagMatrix(CFNULL)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

VolTermDiagBlockJacobSpectralFD::~VolTermDiagBlockJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void VolTermDiagBlockJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void VolTermDiagBlockJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  VolTermRHSJacobSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void VolTermDiagBlockJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  // get the diagonal block Jacobian matrix datahandle
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

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

        // COMPUTE TOTAL RESIDUAL
        m_resUpdates += m_diffResUpdates;

        // COMPUTE JACOBIAN CONTRIBUTION OF VOLUME TERMS
        computeCellGradsMinusBndFaceTerms();
        m_currDiagMatrix = &diagBlockJacobMatr[elemIdx];
        computeJacobVolTerm();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  CFTRACEEND;
 }

//////////////////////////////////////////////////////////////////////////////

void VolTermDiagBlockJacobSpectralFD::setVolumeAndFaceTermComputersData()
{
  // call parent class function
  VolTermRHSSpectralFD::setVolumeAndFaceTermComputersData();
}

//////////////////////////////////////////////////////////////////////////////

void VolTermDiagBlockJacobSpectralFD::addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx)
{
  const CFuint nbrRes = m_currDiagMatrix->nbRows();
  cf_assert(nbrRes == m_derivResUpdates.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    (*m_currDiagMatrix)(iRes,pertResUpdIdx) += m_derivResUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermDiagBlockJacobSpectralFD::computeJacobVolTerm()
{
  // get number of solution points in this cell
  const CFuint nbrSolPnts = m_cellStates->size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // loop over the states in this cell to perturb the states
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

      // add the derivative of the residual updates to diagonal block matrix
      addToDiagBlockJacobMatrix(pertResUpdIdx);

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the volume and face flux points
      restoreVolumeAndFacePhysVars(iEqPert);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolTermDiagBlockJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  VolTermRHSSpectralFD::setup();

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // resize m_physVarGradUpdates
  m_physVarGradUpdates.resize(2);
}

//////////////////////////////////////////////////////////////////////////////

void VolTermDiagBlockJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of the parent class
  VolTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    VolTermDiagBlockJacobSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = VolTermRHSJacobSpectralFD::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr     );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
