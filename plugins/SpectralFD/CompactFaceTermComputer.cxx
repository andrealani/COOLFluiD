#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/CompactFaceTermComputer.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

CompactFaceTermComputer::CompactFaceTermComputer(const std::string& name) :
  BaseFaceTermComputer(name),
  m_faceCompactDiffFluxComputer(CFNULL),
  m_coefSolPolyDerivInFlxPnts(CFNULL),
  m_solPntIdxsSolPolyDerivInFlxPnts(CFNULL),
  m_cellInvJacobMatr(),
  m_auxMatrix(),
  m_matrInverter(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

CompactFaceTermComputer::~CompactFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::configure ( Config::ConfigArgs& args )
{
  BaseFaceTermComputer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::setFaceTermData()
{
  // call setFaceTermData() of parent class
  BaseFaceTermComputer::setFaceTermData();

  // get the SpectralFDElementData
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get coefficients for solution polynomial derivation
  m_coefSolPolyDerivInFlxPnts = sdLocalData[0]->getCoefSolPolyDerivInFlxPntsOptim();

  // get solution point indexes for solution polynomial derivation
  m_solPntIdxsSolPolyDerivInFlxPnts = sdLocalData[0]->getSolPntIdxsSolPolyDerivInFlxPntsOptim();
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::computeNeighbourCellData()
{
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get neighbouring cell
    GeometricEntity* cell = m_face->getNeighborGeo(iSide);

    // compute volume
    m_cellVolumes[iSide] = cell->computeVolume();

    // compute inverse Jacobian matrices and flux points characteristic length scales
    m_cellInvJacobMatr[iSide] = cell->
        computeGeometricShapeFunctionJacobianMatrix
        ((*m_faceFlxPntCellMappedCoords)[m_orient][iSide]);
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      // get Jacobian matrix
      m_auxMatrix = m_cellInvJacobMatr[iSide][iFlx];

      // compute inverse Jacobian matrix
      m_matrInverter->invert(m_auxMatrix,m_cellInvJacobMatr[iSide][iFlx]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::reconstructFluxPntsSolPolyGrads(const vector< vector< State* >* >& cellStates)
{
  cf_assert(cellStates.size() == 2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_statesReconstr->computePolyGradients(*cellStates[iSide],m_flxPntGrads[iSide],
                                           *m_coefSolPolyDerivInFlxPnts,
                                           (*m_faceFlxPntConn)[m_orient][iSide],
                                           *m_solPntIdxsSolPolyDerivInFlxPnts,
                                           m_cellInvJacobMatr[iSide]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::backupAndReconstrFluxPntsSolPolyGrad(const CFuint side,
                                                                   const CFuint iVar,
                                                                   const vector< State* >& cellStates)
{
  // backup
  backupPhysVarGrad(side,iVar);

  // reconstruct
  m_statesReconstr->computePolyGradient(iVar,cellStates,m_flxPntGrads[side],
                                        *m_coefSolPolyDerivInFlxPnts,
                                        (*m_faceFlxPntConn)[m_orient][side],
                                        *m_solPntIdxsSolPolyDerivInFlxPnts,
                                        m_cellInvJacobMatr[side]);
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::computePhysVarGradFaceTerm(const CFuint iVar, vector< vector< RealVector > >& physVarGradUpdates)
{
  // compute the averaged physical variable in the flux points (stored in m_flxPntRiemannFlux)
  // (m_leftRVSol and m_rightRVSol are RealVector pointers to the left and right states)
  m_flxPntRiemannFlux = m_faceCompactDiffFluxComputer->computeAvgGradVar(iVar,m_flxPntRVSol[LEFT],m_flxPntRVSol[RIGHT],
                                                                         m_nbrFlxPnts[m_orient]);

  // compute the actual volume term contribution to the gradient
  computePhysVarGradFaceTermFromFlxPntSol(iVar,physVarGradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::computePhysVarGradFaceTermFromFlxPntSol(const CFuint iVar,
                                                                      vector< vector< RealVector > >& physVarGradUpdates)
{
  cf_assert(physVarGradUpdates.size() == 2);

  // number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // loop over sides
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // set updates to zero
    const CFuint nbrSolPnts = physVarGradUpdates[iSide].size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      physVarGradUpdates[iSide][iSol] = 0.0;
    }

    // loop over flux points
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      // flux point index
      const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iSide][iFlx];

      // flux point index in the matrix m_solPntsDerivCoefs
      const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

      // compute gradient term
      m_gradTerm = m_flxPntRiemannFlux[iFlx][iVar]*
                   m_faceJacobVecSizeFlxPnts[iFlx][iSide]*
                   m_unitNormalFlxPnts[iFlx];

      // add contribution of this flux point to the solution points
      for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
      {
        // solution point index
        const CFuint solIdx = (*m_solPntIdxsForDerivation)[flxIdx][iSol];
        physVarGradUpdates[iSide][solIdx] +=
            (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*m_gradTerm;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  BaseFaceTermComputer::setup();

  // cast face diffusive flux computer to IP face diffusive flux computer
  m_faceCompactDiffFluxComputer = m_faceDiffFluxComputer.d_castTo<FaceDiffusiveFluxCompactApproach>();

  // create matrix inverter
  m_matrInverter = MathTools::MatrixInverter::create(m_dim,false);

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of solution points
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize variables
  m_cellInvJacobMatr.resize(2);
  m_cellInvJacobMatr[LEFT ].resize(nbrSolPnts,RealMatrix(m_dim,m_dim));
  m_cellInvJacobMatr[RIGHT].resize(nbrSolPnts,RealMatrix(m_dim,m_dim));
  m_auxMatrix.resize(m_dim,m_dim);
}

//////////////////////////////////////////////////////////////////////////////

void CompactFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_matrInverter);

  // call unsetup of parent class
  BaseFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
