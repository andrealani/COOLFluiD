#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxBR2ApproachLES.hh"
#include "SpectralFDLES/LESBR2FaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LESBR2FaceTermComputer,SpectralFDMethodData,BaseFaceTermComputer,SpectralFDLESModule >
LESBR2FaceTermComputerProvider("LESBR2FaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

LESBR2FaceTermComputer::LESBR2FaceTermComputer(const std::string& name) :
  NSBR2FaceTermComputer(name),
  m_filterWidthVolumes()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LESBR2FaceTermComputer::~LESBR2FaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void LESBR2FaceTermComputer::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NSBR2FaceTermComputer::setup();

}

//////////////////////////////////////////////////////////////////////////////

void LESBR2FaceTermComputer::setFaceTermData()
{
  // Call setFaceTermData() of the parent (but which is defined in BR2FaceTermComputer)
  NSBR2FaceTermComputer::setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void LESBR2FaceTermComputer::computeNeighbourCellData()
{
  // Call computeNeighbourCellData() of parent class
  NSBR2FaceTermComputer::computeNeighbourCellData();

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

void LESBR2FaceTermComputer::computeDiffFaceTerm(vector< RealVector>& resUpdates)
{
  CFAUTOTRACE;

  // Set the filter width volume in the facediffusiveflux computer
  m_faceDiffFluxComputer.d_castTo<FaceDiffusiveFluxBR2ApproachLES>()->setFilterWidthVolumes(&m_filterWidthVolumes);

  // Call computeDiffFaceTerm(resUpdates) of the parent (but which is defined in BaseFaceTermComputer)
  NSBR2FaceTermComputer::computeDiffFaceTerm(resUpdates);

}

//////////////////////////////////////////////////////////////////////////////

void LESBR2FaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NSBR2FaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
