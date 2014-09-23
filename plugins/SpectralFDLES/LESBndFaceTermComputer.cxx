#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/NavierStokesBndFaceTermComputer.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxLocalApproachLES.hh"
#include "SpectralFDLES/LESBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LESBndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDLESModule >
LESBndFaceTermComputerProvider("LESBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

LESBndFaceTermComputer::LESBndFaceTermComputer(const std::string& name) : NavierStokesBndFaceTermComputer(name),
  m_filterWidthVolumes()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LESBndFaceTermComputer::~LESBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void LESBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NavierStokesBndFaceTermComputer::setup();
}

//////////////////////////////////////////////////////////////////////////////

void LESBndFaceTermComputer::setFaceTermData()
{
  CFAUTOTRACE;

  // Call setFaceTermData() of the parent class
  NavierStokesBndFaceTermComputer::setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void LESBndFaceTermComputer::computeNeighbourCellData()
{
  CFAUTOTRACE;

  // Call computeNeighbourCellData() of parent class
  NavierStokesBndFaceTermComputer::computeNeighbourCellData();

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

void LESBndFaceTermComputer::computeDiffFaceTerm(RealVector& resUpdates)
{
  CFAUTOTRACE;

  // Set the filter width volume in the facediffusiveflux computer
  m_faceDiffFluxComputer.d_castTo<FaceDiffusiveFluxLocalApproachLES>()->setFilterWidthVolumes(&m_filterWidthVolumes);

  // Call computeDiffFaceTerm(resUpdates) of the parent class
  NavierStokesBndFaceTermComputer::computeDiffFaceTerm(resUpdates);

}

//////////////////////////////////////////////////////////////////////////////

void LESBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NavierStokesBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
