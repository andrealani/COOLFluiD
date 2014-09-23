#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "SpectralFD/FaceTermDiagBlockJacobSpectralFD.hh"
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

MethodCommandProvider<FaceTermDiagBlockJacobSpectralFD, SpectralFDMethodData, SpectralFDModule>
  FaceTermDiagBlockJacobSpectralFDProvider("FaceTermDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

FaceTermDiagBlockJacobSpectralFD::FaceTermDiagBlockJacobSpectralFD(const std::string& name) :
  FaceTermRHSJacobSpectralFD(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  m_currDiagMatrix()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

FaceTermDiagBlockJacobSpectralFD::~FaceTermDiagBlockJacobSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermDiagBlockJacobSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermDiagBlockJacobSpectralFD::configure ( Config::ConfigArgs& args )
{
  FaceTermRHSJacobSpectralFD::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermDiagBlockJacobSpectralFD::execute()
{
  CFTRACEBEGIN;

  /// @note in this function it is assumed that there is only one type of cell neighbouring the faces

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

  // SET THE DATA NEEDED TO COMPUTE THE FACE TERMS
  setFaceTermComputerData();

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // SET THE ORIENTATION OF THE FACES
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

      // if one of the neighbouring cells is parallel updatable, compute the face terms
      const bool lParUpdatable = (*m_states[LEFT ])[0]->isParUpdatable();
      const bool rParUpdatable = (*m_states[RIGHT])[0]->isParUpdatable();
      if (lParUpdatable || rParUpdatable)
      {
        // COMPUTE DATA IN FACE TERM COMPUTER
        computeFaceTermData();

        // COMPUTE CONVECTIVE FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
        m_faceTermComputer->
            computeConvFaceTermAndWaveSpeedUpdates(m_resUpdates,m_updateCoefContr);

        // COMPUTE DIFFUSIVE FACE TERMS AND UPDATE COEFFICIENT CONTRIBUTIONS
        m_faceTermComputer->
            computeDiffFaceTermAndUpdateCoefContributions(m_diffResUpdates,m_diffUpdateCoefContr);

        // ADD TOTAL UPDATE TO UPDATE COEFFICIENTS
        m_resUpdates[LEFT ] += m_diffResUpdates[LEFT ];
        m_resUpdates[RIGHT] += m_diffResUpdates[RIGHT];
        m_updateCoefContr[LEFT ] += m_diffUpdateCoefContr[LEFT ];
        m_updateCoefContr[RIGHT] += m_diffUpdateCoefContr[RIGHT];
        addUpdateCoeffContributions();

        // COMPUTE JACOBIAN CONTRIBUTION OF FACE TERMS
        if (lParUpdatable)
        {
          m_currDiagMatrix = &diagBlockJacobMatr[m_cells[LEFT ]->getID()];
          computeOneJacobFaceTerm(LEFT);
        }
        if (rParUpdatable)
        {
          m_currDiagMatrix = &diagBlockJacobMatr[m_cells[RIGHT]->getID()];
          computeOneJacobFaceTerm(RIGHT);
        }
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
// CF_DEBUG_EXIT;
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermDiagBlockJacobSpectralFD::computeOneJacobFaceTerm(const CFuint side)
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

      // backup and reconstruct physical variable and its gradient in the flux points
      backupAndReconstructFacePhysVarsAndGrad(side,iEqPert);

      // compute the perturbed face term
      m_faceTermComputer->computeConvFaceTerm(m_pertResUpdates);
      m_faceTermComputer->computeDiffFaceTerm(m_diffResUpdates);

      // add contributions to the Jacobian
      // compute the finite difference derivative of the face term
      m_pertResUpdates[side] += m_diffResUpdates[side];
      m_numJacob->computeDerivative(m_pertResUpdates[side],m_resUpdates[side],m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to diagonal block matrix
      addToDiagBlockJacobMatrix(pertResUpdIdx);

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable and its gradient in the flux points
      restoreFacePhysVarsAndGrad(side,iEqPert);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermDiagBlockJacobSpectralFD::addToDiagBlockJacobMatrix(const CFuint pertResUpdIdx)
{
  const CFuint nbrRes = m_currDiagMatrix->nbRows();
  cf_assert(nbrRes <= m_derivResUpdates.size());
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes)
  {
    (*m_currDiagMatrix)(iRes,pertResUpdIdx) += m_derivResUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void FaceTermDiagBlockJacobSpectralFD::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  FaceTermRHSSpectralFD::setup();

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

void FaceTermDiagBlockJacobSpectralFD::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of parent class
  FaceTermRHSSpectralFD::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    FaceTermDiagBlockJacobSpectralFD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = FaceTermRHSJacobSpectralFD::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
