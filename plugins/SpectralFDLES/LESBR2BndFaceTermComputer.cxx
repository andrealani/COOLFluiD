#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/NSBR2BndFaceTermComputer.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxBR2ApproachLES.hh"
#include "SpectralFDLES/LESBR2BndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LESBR2BndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDLESModule >
LESBR2BndFaceTermComputerProvider("LESBR2BndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

LESBR2BndFaceTermComputer::LESBR2BndFaceTermComputer(const std::string& name) : NSBR2BndFaceTermComputer(name),
  m_filterWidthVolumes()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LESBR2BndFaceTermComputer::~LESBR2BndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void LESBR2BndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NSBR2BndFaceTermComputer::setup();
}

//////////////////////////////////////////////////////////////////////////////

void LESBR2BndFaceTermComputer::setFaceTermData()
{
  CFAUTOTRACE;

  // Call setFaceTermData() of the parent (but which is defined in BR2BndFaceTermComputer)
  NSBR2BndFaceTermComputer::setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void LESBR2BndFaceTermComputer::computeNeighbourCellData()
{
  CFAUTOTRACE;

  // Call computeNeighbourCellData() of parent class
  NSBR2BndFaceTermComputer::computeNeighbourCellData();

  // Resize m_filterWidthVolumes
  m_filterWidthVolumes.resize(m_nbrFlxPnts[m_orient]);

  // Get neighbouring cell
  GeometricEntity* cell = m_face->getNeighborGeo(0);

  // Compute Jacobian determinants
  std::valarray<CFreal> jacobDets(m_nbrFlxPnts[m_orient]);
  jacobDets = cell->computeGeometricShapeFunctionJacobianDeterminant((*m_faceFlxPntCellMappedCoords)[m_orient]);

  // Compute inverse characteristic lengths
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_filterWidthVolumes[iFlx] = jacobDets[iFlx];
  }

}

//////////////////////////////////////////////////////////////////////////////

void LESBR2BndFaceTermComputer::computeDiffFaceTerm(RealVector& resUpdates)
{
  CFAUTOTRACE;

  // Set the filter width volume in the facediffusiveflux computer
  m_faceDiffFluxComputer.d_castTo<FaceDiffusiveFluxBR2ApproachLES>()->setFilterWidthVolumes(&m_filterWidthVolumes);

  // Call computeDiffFaceTerm(resUpdates) of the parent class (but which is defined in BaseFaceTermComputer)
  NSBR2BndFaceTermComputer::computeDiffFaceTerm(resUpdates);

}

//////////////////////////////////////////////////////////////////////////////

void LESBR2BndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NSBR2BndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
