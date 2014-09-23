#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/CompactBndFaceTermComputer.hh"
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

CompactBndFaceTermComputer::CompactBndFaceTermComputer(const std::string& name) :
  BaseBndFaceTermComputer(name),
  m_coefSolPolyDerivInFlxPnts(CFNULL),
  m_solPntIdxsSolPolyDerivInFlxPnts(CFNULL),
  m_cellInvJacobMatr(),
  m_auxMatrix(),
  m_matrInverter(),
  m_solPolyDerivCoefsPntSet(),
  m_cellInvJacobMatrPntSet()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

CompactBndFaceTermComputer::~CompactBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::configure ( Config::ConfigArgs& args )
{
  BaseBndFaceTermComputer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::setFaceTermData()
{
  // call setFaceTermData of parent class
  BaseBndFaceTermComputer::setFaceTermData();

  // get the SpectralFDElementData
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get coefficients for solution polynomial derivation
  m_coefSolPolyDerivInFlxPnts = sdLocalData[0]->getCoefSolPolyDerivInFlxPntsOptim();

  // get solution point indexes for solution polynomial derivation
  m_solPntIdxsSolPolyDerivInFlxPnts = sdLocalData[0]->getSolPntIdxsSolPolyDerivInFlxPntsOptim();
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::computeNeighbourCellData()
{
  // get neighbouring cell
  GeometricEntity* cell = m_face->getNeighborGeo(0);

  // compute volume
  m_cellVolume = cell->computeVolume();

  // compute inverse Jacobian matrices and flux points characteristic length scales
  m_cellInvJacobMatr = cell->
      computeGeometricShapeFunctionJacobianMatrix((*m_faceFlxPntCellMappedCoords)[m_orient]);
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    // get Jacobian matrix
    m_auxMatrix = m_cellInvJacobMatr[iFlx];

    // compute inverse Jacobian matrix
    m_matrInverter->invert(m_auxMatrix,m_cellInvJacobMatr[iFlx]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::reconstructFluxPntsSolPolyGrads(const vector< State* >& cellIntStates)
{
  // compute the internal gradients
  m_statesReconstr->computePolyGradients(cellIntStates,m_flxPntIntGrads,
                                         *m_coefSolPolyDerivInFlxPnts,
                                         (*m_faceFlxPntConn)[m_orient],
                                         *m_solPntIdxsSolPolyDerivInFlxPnts,
                                         m_cellInvJacobMatr);

  // compute the ghost gradients
  computeGhostGradients();
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::backupAndReconstrFluxPntsSolPolyGrad(const CFuint iVar,
                                                                      const vector< State* >& cellIntStates)
{
  // backup
  backupGradPhysVar(iVar);

  // reconstruct
  m_statesReconstr->computePolyGradient(iVar,cellIntStates,m_flxPntIntGrads,
                                        *m_coefSolPolyDerivInFlxPnts,
                                        (*m_faceFlxPntConn)[m_orient],
                                        *m_solPntIdxsSolPolyDerivInFlxPnts,
                                        m_cellInvJacobMatr);

  // compute the ghost gradients
  computeGhostGradients();
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::setPointSet(const std::vector< RealVector >& faceMappedCoordPntSet)
{
  // call setPointSet of parent class
  BaseBndFaceTermComputer::setPointSet(faceMappedCoordPntSet);

  // get the SpectralFDElementData
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // compute solution polynomial derivation coefficients in given points set
  m_solPolyDerivCoefsPntSet = sdLocalData[0]->getSolPolyDerivsAtNode(m_cellMappedCoordsPntSet);

  // number of points in point set
  const CFuint nbrPnts = faceMappedCoordPntSet.size();

  // resize some variables
  if (nbrPnts > m_cellInvJacobMatrPntSet.size())
  {
    m_cellInvJacobMatrPntSet.resize(nbrPnts,RealMatrix(m_dim,m_dim));
  }
  else if (nbrPnts < m_cellInvJacobMatrPntSet.size())
  {
    m_cellInvJacobMatrPntSet.resize(nbrPnts,RealMatrix(m_dim,m_dim));
  }
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::computeFacePntSetData()
{
  // call computeFaceOutputData() of parent class
  BaseBndFaceTermComputer::computeFacePntSetData();

  // get neighbouring cell
  GeometricEntity* cell = m_face->getNeighborGeo(0);

  // compute inverse Jacobian matrices and flux points characteristic length scales
  m_cellInvJacobMatrPntSet =
      cell->computeGeometricShapeFunctionJacobianMatrix(m_cellMappedCoordsPntSet);
  const CFuint nbrPnts = m_cellMappedCoordsPntSet.size();
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    // get Jacobian matrix
    m_auxMatrix = m_cellInvJacobMatrPntSet[iPnt];

    // compute inverse Jacobian matrix
    m_matrInverter->invert(m_auxMatrix,m_cellInvJacobMatrPntSet[iPnt]);
  }
}

//////////////////////////////////////////////////////////////////////////////

vector< vector< RealVector > >& CompactBndFaceTermComputer::reconstructGivenPntsGrads(const vector< State* >& cellIntStates)
{
  // number of points
  const CFuint nbrPnts = m_solPolyDerivCoefsPntSet.size();
  cf_assert(nbrPnts == m_gradsPntSet.size());

  // number of solution points/basis polynomials
  const CFuint nbrPolys = cellIntStates.size();
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      RealVector grad(m_dim);
      grad = 0.0;
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        cf_assert(nbrPolys == m_solPolyDerivCoefsPntSet[iPnt][iDim].size());
        for (CFuint iPoly = 0; iPoly < nbrPolys; ++iPoly)
        {
          grad[iDim] += m_solPolyDerivCoefsPntSet[iPnt][iDim][iPoly]*(*cellIntStates[iPoly])[iEq];
        }
      }

      // transform the gradient
      *m_intGradsPntSet[iPnt][iEq] = m_cellInvJacobMatrPntSet[iPnt]*grad;
    }
  }

  // compute ghost gradients
  m_bcStateComputer->computeGhostGradients(m_intGradsPntSet,m_ghostGradsPntSet,m_unitNormalPntSet,m_coordsPntSet);

  // average internal and ghost gradients
  for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsPntSet[iPnt][iEq] = 0.5*(*m_intGradsPntSet[iPnt][iEq] + *m_ghostGradsPntSet[iPnt][iEq]);
    }
  }

  /// @todo: add lifting operators

  return m_gradsPntSet;
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  BaseBndFaceTermComputer::setup();

  // create matrix inverter
  m_matrInverter = MathTools::MatrixInverter::create(m_dim,false);

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of solution points
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize m_cellInvJacobMatr
  m_cellInvJacobMatr.resize(nbrSolPnts,RealMatrix(m_dim,m_dim));
  m_auxMatrix.resize(m_dim,m_dim);
}

//////////////////////////////////////////////////////////////////////////////

void CompactBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_matrInverter);

  // call unsetup of parent class
  BaseBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
