#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/IPFaceTermComputer.hh"
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
  IPFaceTermComputer,SpectralFDMethodData,BaseFaceTermComputer,SpectralFDModule >
  IPFaceTermComputerProvider("IPFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

IPFaceTermComputer::IPFaceTermComputer(const std::string& name) :
  CompactFaceTermComputer(name),
  m_backupPhysVarGrad()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

IPFaceTermComputer::~IPFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void IPFaceTermComputer::configure ( Config::ConfigArgs& args )
{
  CompactFaceTermComputer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void IPFaceTermComputer::computeNeighbourCellData()
{
  // set face characteristic lengths to zero
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_faceInvCharLengths[iFlx] = 0.0;
  }

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

      // contribution to characteristic lengths
      m_faceInvCharLengths[iFlx] += m_dim == 2 ? m_auxMatrix.determ2() : m_auxMatrix.determ3();

      // compute inverse Jacobian matrix
      m_matrInverter->invert(m_auxMatrix,m_cellInvJacobMatr[iSide][iFlx]);
    }
  }

  // compute inverse characteristic lengths
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_faceInvCharLengths[iFlx] = m_faceJacobVecAbsSizeFlxPnts[iFlx]/m_faceInvCharLengths[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void IPFaceTermComputer::backupPhysVarGrad(const CFuint side, const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVarGrad.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_backupPhysVarGrad[iFlx] = *m_flxPntGrads[side][iFlx][iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void IPFaceTermComputer::restorePhysVarGrad(const CFuint side, const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVarGrad.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    *m_flxPntGrads[side][iFlx][iVar] = m_backupPhysVarGrad[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void IPFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  CompactFaceTermComputer::setup();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of flux points
  const CFuint nbrFlxPnts = sdLocalData[0]->getNbrOfFaceFlxPnts();

  // resize m_backupPhysVarGrad
  m_backupPhysVarGrad.resize(nbrFlxPnts,RealVector(m_dim));
}

//////////////////////////////////////////////////////////////////////////////

void IPFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // call unsetup of parent class
  CompactFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
