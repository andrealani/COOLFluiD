#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxIPApproachLES.hh"
#include "SpectralFDLES/LESIPFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LESIPFaceTermComputer,SpectralFDMethodData,BaseFaceTermComputer,SpectralFDLESModule >
LESIPFaceTermComputerProvider("LESIPFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

LESIPFaceTermComputer::LESIPFaceTermComputer(const std::string& name) :
  NSIPFaceTermComputer(name),
  m_filterWidthVolumes()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LESIPFaceTermComputer::~LESIPFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void LESIPFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NSIPFaceTermComputer::setup();

}

//////////////////////////////////////////////////////////////////////////////

void LESIPFaceTermComputer::setFaceTermData()
{
  // Call setFaceTermData() of the parent (but which is defined in IPFaceTermComputer)
  NSIPFaceTermComputer::setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void LESIPFaceTermComputer::computeNeighbourCellData()
{
  // Call computeNeighbourCellData() of parent class
  NSIPFaceTermComputer::computeNeighbourCellData();

  // Resize m_filterWidthVolumes
  m_filterWidthVolumes.resize(m_nbrFlxPnts[m_orient]);

  vector< std::valarray<CFreal> > jacobDets(2,std::valarray<CFreal>(m_nbrFlxPnts[m_orient]));
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // Get neighbouring cell
    GeometricEntity* cell = m_face->getNeighborGeo(iSide);

    // Compute Jacobian determinants
    jacobDets[iSide] = cell->computeGeometricShapeFunctionJacobianDeterminant((*m_faceFlxPntCellMappedCoords)[m_orient][iSide]);
  }

  // Set averaged jacobian determinants in the LESVarSet
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    // Compute m_filterWidthVolumes
    const CFreal avgJacobDets = 0.5*(jacobDets[LEFT][iFlx] + jacobDets[RIGHT][iFlx]);

    m_filterWidthVolumes[iFlx] = avgJacobDets;
  }


}

//////////////////////////////////////////////////////////////////////////////

void LESIPFaceTermComputer::computeDiffFaceTerm(vector< RealVector>& resUpdates)
{
  CFAUTOTRACE;

  // Set the filter width volume in the facediffusiveflux computer
  m_faceDiffFluxComputer.d_castTo<FaceDiffusiveFluxIPApproachLES>()->setFilterWidthVolumes(&m_filterWidthVolumes);

  // Call computeDiffFaceTerm(resUpdates) of the parent (but which is defined in BaseFaceTermComputer)
  NSIPFaceTermComputer::computeDiffFaceTerm(resUpdates);

}

//////////////////////////////////////////////////////////////////////////////

void LESIPFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NSIPFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
