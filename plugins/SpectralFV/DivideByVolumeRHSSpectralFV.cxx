#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/DivideByVolumeRHSSpectralFV.hh"
#include "SpectralFV/SpectralFVElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DivideByVolumeRHSSpectralFV, SpectralFVMethodData, SpectralFVModule> DivideByVolumeRHSSpectralFVProvider("DivideByVolumeRHS");

//////////////////////////////////////////////////////////////////////////////

DivideByVolumeRHSSpectralFV::DivideByVolumeRHSSpectralFV(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  m_cellBuilder(CFNULL),
  m_invVolFracCVs(),
  m_cell(),
  m_cellStates(),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

DivideByVolumeRHSSpectralFV::~DivideByVolumeRHSSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void DivideByVolumeRHSSpectralFV::configure ( Config::ConfigArgs& args )
{
  SpectralFVMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DivideByVolumeRHSSpectralFV::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // get the local spectral FV data
    vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

    // get the inverse volume fractions of the CVs
    m_invVolFracCVs = svLocalData[iElemType]->getInvVolFracCV();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // divide residuals by cell volume
      divideRHSByVolume();

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void DivideByVolumeRHSSpectralFV::divideRHSByVolume()
{
  // get the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get the updateCoeff
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();

  // get the cell volume
  const CFreal invCellVolume = 1.0/m_cell->computeVolume();

  // get number of CVs
  const CFuint nbrCVs = m_cellStates->size();

  // loop over residuals
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    // get CV ID
    const CFuint cvID = (*m_cellStates)[iCV]->getLocalID();

    // inverse CV volume
    const CFreal invVolume = invCellVolume*(*m_invVolFracCVs)[iCV];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // divide rhs by volume
      rhs(cvID,iEq,m_nbrEqs) *= invVolume;
    }

    // divide update coeff by volume
    updateCoeff[cvID] *= invVolume;
  }
}

//////////////////////////////////////////////////////////////////////////////

void DivideByVolumeRHSSpectralFV::setup()
{
  CFAUTOTRACE;

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();
}

//////////////////////////////////////////////////////////////////////////////

void DivideByVolumeRHSSpectralFV::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
DivideByVolumeRHSSpectralFV::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
