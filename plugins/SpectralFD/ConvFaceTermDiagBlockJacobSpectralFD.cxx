#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/ConvFaceTermDiagBlockJacobSpectralFD.hh"
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

MethodCommandProvider<ConvFaceTermDiagBlockJacobSpectralFD, SpectralFDMethodData, SpectralFDModule> ConvFaceTermDiagBlockJacobSpectralFDProvider("ConvFaceTermDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermDiagBlockJacobSpectralFD::ConvFaceTermDiagBlockJacobSpectralFD(const std::string& name) :
    ConvFaceTermRHSJacobSpectralFD(name),
    socket_diagBlockJacobMatr("diagBlockJacobMatr"),
    m_currDiagMatrix()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ConvFaceTermDiagBlockJacobSpectralFD::~ConvFaceTermDiagBlockJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermDiagBlockJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermDiagBlockJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  ConvFaceTermRHSJacobSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermDiagBlockJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  /// @note in this function it is assumed that there is only one type of cell neighbouring the faces

  // set the data needed to compute the face terms;
  setFaceTermData();

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

  // get the diagonal block Jacobian matrix datahandle
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = faces;
  geoData.isBoundary = false;

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // set the orientation of the faces
    m_faceTermComputer->setFaceOrientation(m_orient);

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoData.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      // if one of the neighbouring cells is parallel updatable or the gradients have to be computed,
      // set face data and reconstruct states
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable() || hasDiffTerm)
      {
        // set the current face and compute the face data in the face term computer
        m_faceTermComputer->setCurrentFace(m_face);
        m_faceTermComputer->computeFaceData();

        // reconstruct the states
        m_faceTermComputer->reconstructFluxPntsStates(m_states);
      }

      // if one of the neighbouring cells is parallel updatable,
      // compute the face term
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
        // compute the face terms and the wave speed updates
        m_faceTermComputer
            ->computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_waveSpeedUpd);

        // update the wave speeds
        updateWaveSpeed();

        if ((*m_states[LEFT ])[0]->isParUpdatable())
        {
          m_currDiagMatrix = &diagBlockJacobMatr[m_cells[LEFT ]->getID()];
          computeOneDiagBlockJacobConvFaceTerm(LEFT );
        }
        if ((*m_states[RIGHT])[0]->isParUpdatable())
        {
          m_currDiagMatrix = &diagBlockJacobMatr[m_cells[RIGHT]->getID()];
          computeOneDiagBlockJacobConvFaceTerm(RIGHT);
        }
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradientFaceTerm();
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermDiagBlockJacobSpectralFD::computeOneDiagBlockJacobConvFaceTerm(const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // loop over the states to perturb the states
  const CFuint nbrSolPnts = m_states[side]->size();
  CFuint pertResUpdIdx = 0;
  for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
  {
    // dereference states
    State& pertState = *(*m_states[side])[iSolPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert, ++pertResUpdIdx)
    {
      // perturb physical variable in state
      m_numJacob->perturb(iEqPert,pertState[iEqPert]);

      // backup and reconstruct physical variable in the flux points
      m_faceTermComputer->backupAndReconstructPhysVar(side,iEqPert,*m_states[side]);

      // compute the perturbed face term
      m_faceTermComputer->computeConvFaceTerm(m_pertResUpdates);

      // add contributions to the Jacobian
      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_pertResUpdates[side],m_resUpdates[side],m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to diagonal block matrix
      addToDiagBlockJacobMatrix(pertResUpdIdx);

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the flux points
      m_faceTermComputer->restorePhysVar(side,iEqPert);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermDiagBlockJacobSpectralFD::addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx)
{
  const CFuint nbrRes = m_currDiagMatrix->nbRows();
  cf_assert(nbrRes <= m_derivResUpdates.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    (*m_currDiagMatrix)(iRes,pertResUpdIdx) += m_derivResUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermDiagBlockJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // call setup of parent class of parent class, in order not to get the lss
  ConvFaceTermRHSSpectralFD::setup();

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // maximum number of solution points in a cell
  const CFuint maxNbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize variables
  const CFuint resSize = maxNbrSolPnts*m_nbrEqs;
  m_pertResUpdates.resize(2);
  m_pertResUpdates[LEFT ].resize(resSize);
  m_pertResUpdates[RIGHT].resize(resSize);
  m_derivResUpdates.resize(resSize);
}

//////////////////////////////////////////////////////////////////////////////

void ConvFaceTermDiagBlockJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of parent class of parent class, in order not to get the lss
  ConvFaceTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    ConvFaceTermDiagBlockJacobSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ConvFaceTermRHSJacobSpectralFD::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
