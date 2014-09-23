#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/IPBndFaceTermComputer.hh"
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

Framework::MethodStrategyProvider<
    IPBndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDModule >
IPBndFaceTermComputerProvider("IPBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

IPBndFaceTermComputer::IPBndFaceTermComputer(const std::string& name) :
  CompactBndFaceTermComputer(name),
  m_backupPhysVarGrad(),
  m_backupGhostStateGrad()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

IPBndFaceTermComputer::~IPBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void IPBndFaceTermComputer::configure ( Config::ConfigArgs& args )
{
  CompactBndFaceTermComputer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void IPBndFaceTermComputer::computeNeighbourCellData()
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

    // compute inverse characteristic length scale
    m_faceInvCharLengths[iFlx] =
        0.5*m_faceJacobVecAbsSizeFlxPnts[iFlx]/(m_dim == 2 ? m_auxMatrix.determ2() : m_auxMatrix.determ3());

    // compute inverse Jacobian matrix
    m_matrInverter->invert(m_auxMatrix,m_cellInvJacobMatr[iFlx]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void IPBndFaceTermComputer::backupGradPhysVar(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVarGrad.size());
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStateGrad.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_backupPhysVarGrad[iFlx] = *m_flxPntIntGrads[iFlx][iVar];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_backupGhostStateGrad[iFlx][iEq] = *m_flxPntGhostGrads[iFlx][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void IPBndFaceTermComputer::restorePhysVarGrad(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVarGrad.size());
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStateGrad.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    *m_flxPntIntGrads[iFlx][iVar] = m_backupPhysVarGrad[iFlx];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *m_flxPntGhostGrads[iFlx][iEq] = m_backupGhostStateGrad[iFlx][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void IPBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  CompactBndFaceTermComputer::setup();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of flux points
  const CFuint nbrFlxPnts = sdLocalData[0]->getNbrOfFaceFlxPnts();

  // resize m_backupPhysVarGrad and m_backupGhostStateGrad
  m_backupPhysVarGrad.resize(nbrFlxPnts,RealVector(m_dim));
  m_backupGhostStateGrad.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_backupGhostStateGrad[iFlx].resize(m_nbrEqs,RealVector(m_dim));
  }
}

//////////////////////////////////////////////////////////////////////////////

void IPBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of parent class
  CompactBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
