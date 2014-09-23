#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/NSIPBndFaceTermComputer.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxIPApproachLES.hh"
#include "SpectralFDLES/LESIPBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LESIPBndFaceTermComputer,SpectralFDMethodData,BaseBndFaceTermComputer,SpectralFDLESModule >
LESIPBndFaceTermComputerProvider("LESIPBndFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

LESIPBndFaceTermComputer::LESIPBndFaceTermComputer(const std::string& name) :
  NSIPBndFaceTermComputer(name),
  m_filterWidthVolumes()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LESIPBndFaceTermComputer::~LESIPBndFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void LESIPBndFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NSIPBndFaceTermComputer::setup();
}

//////////////////////////////////////////////////////////////////////////////

void LESIPBndFaceTermComputer::setFaceTermData()
{
  CFAUTOTRACE;

  // Call setFaceTermData() of the parent (but which is defined in IPBndFaceTermComputer)
  NSIPBndFaceTermComputer::setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void LESIPBndFaceTermComputer::computeNeighbourCellData()
{
  CFAUTOTRACE;

  // Call computeNeighbourCellData() of parent class
  NSIPBndFaceTermComputer::computeNeighbourCellData();

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

void LESIPBndFaceTermComputer::computeDiffFaceTerm(RealVector& resUpdates)
{
  CFAUTOTRACE;

  // Set the filter width volume in the facediffusiveflux computer
  m_faceDiffFluxComputer.d_castTo<FaceDiffusiveFluxIPApproachLES>()->setFilterWidthVolumes(&m_filterWidthVolumes);

  // Call computeDiffFaceTerm(resUpdates) of the parent class (but which is defined in BaseFaceTermComputer)
  NSIPBndFaceTermComputer::computeDiffFaceTerm(resUpdates);

}

//////////////////////////////////////////////////////////////////////////////

void LESIPBndFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NSIPBndFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
