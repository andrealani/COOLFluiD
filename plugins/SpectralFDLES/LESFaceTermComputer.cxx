#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/NavierStokesFaceTermComputer.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxLocalApproachLES.hh"
#include "SpectralFDLES/LESFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LESFaceTermComputer,SpectralFDMethodData,BaseFaceTermComputer,SpectralFDLESModule >
LESFaceTermComputerProvider("LESFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

LESFaceTermComputer::LESFaceTermComputer(const std::string& name) :
  NavierStokesFaceTermComputer(name),
  m_filterWidthVolumes()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LESFaceTermComputer::~LESFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void LESFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NavierStokesFaceTermComputer::setup();

}

//////////////////////////////////////////////////////////////////////////////

void LESFaceTermComputer::setFaceTermData()
{
  CFAUTOTRACE;

  // Call setFaceTermData() of the parent class
  NavierStokesFaceTermComputer::setFaceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void LESFaceTermComputer::computeNeighbourCellData()
{
  CFAUTOTRACE;

  // Call computeNeighbourCellData() of parent class
  NavierStokesFaceTermComputer::computeNeighbourCellData();

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
    //cout<< jacobDets[LEFT][iFlx] << endl;

    // Compute m_filterWidthVolumes
    const CFreal avgJacobDets = 0.5*(jacobDets[LEFT][iFlx] + jacobDets[RIGHT][iFlx]);

    m_filterWidthVolumes[iFlx] = avgJacobDets;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LESFaceTermComputer::computeDiffFaceTerm(vector< RealVector>& resUpdates)
{
  CFAUTOTRACE;

  // Set the filter width volume in the facediffusiveflux computer
  m_faceDiffFluxComputer.d_castTo<FaceDiffusiveFluxLocalApproachLES>()->setFilterWidthVolumes(&m_filterWidthVolumes);

  // Call computeDiffFaceTerm(resUpdates) of the parent class
  NavierStokesFaceTermComputer::computeDiffFaceTerm(resUpdates);

}

//////////////////////////////////////////////////////////////////////////////

void LESFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NavierStokesFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
