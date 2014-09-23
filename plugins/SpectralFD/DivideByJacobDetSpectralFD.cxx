#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/DivideByJacobDetSpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DivideByJacobDetSpectralFD, SpectralFDMethodData, SpectralFDModule> DivideByJacobDetSpectralFDProvider("DivideByJacobDet");

//////////////////////////////////////////////////////////////////////////////

DivideByJacobDetSpectralFD::DivideByJacobDetSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_updateCoeff("updateCoeff"),
  m_cellBuilder(CFNULL),
  m_cell(),
  m_solPntsLocalCoords(),
  m_cellStates(),
  m_nbrEqs()
{
}

//////////////////////////////////////////////////////////////////////////////

DivideByJacobDetSpectralFD::~DivideByJacobDetSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void DivideByJacobDetSpectralFD::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DivideByJacobDetSpectralFD::execute()
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

    // get the local spectral FD data
    vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

    // get solution point local coordinates
    m_solPntsLocalCoords = sdLocalData[iElemType]->getSolPntsLocalCoords();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // divide residuals by cell volume
      divideRHSByJacobDet();

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void DivideByJacobDetSpectralFD::divideRHSByJacobDet()
{
  // this is achieved by multiplying the update coefs with the Jacobian determinant
  // (and dividing by the cell volume)

  // get the updateCoeff
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();

  // get the cell volume
  const CFreal invCellVolume = 1.0/m_cell->computeVolume();

  // get jacobian determinants at solution points
  const std::valarray<CFreal> jacobDet =
      m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  // get number of solution points
  const CFuint nbrSolPnts = m_cellStates->size();

  // loop over residuals
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get solution point ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // divide update coeff by volume
    updateCoeff[solID] *= jacobDet[iSol]*invCellVolume;
  }
}

//////////////////////////////////////////////////////////////////////////////

void DivideByJacobDetSpectralFD::setup()
{
  CFAUTOTRACE;

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();
}

//////////////////////////////////////////////////////////////////////////////

void DivideByJacobDetSpectralFD::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
DivideByJacobDetSpectralFD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
