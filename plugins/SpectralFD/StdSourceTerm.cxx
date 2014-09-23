#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"

#include "Framework/NamespaceSwitcher.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/StdSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

// MethodCommandProvider<StdSourceTerm, SpectralFDMethodData, SpectralFDModule>
// StdSourceTermProvider("StdSourceTerm");

//////////////////////////////////////////////////////////////////////////////

StdSourceTerm::StdSourceTerm(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_rhs("rhs"),
  m_cell(CFNULL),
  m_cellStates(CFNULL),
  m_nbrEqs(),
  m_iElemType(),
  m_solPntsLocalCoords(CFNULL),
  m_solPntJacobDets()
{
}

//////////////////////////////////////////////////////////////////////////////

StdSourceTerm::~StdSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::setup()
{
  CFAUTOTRACE;
  SpectralFDMethodCom::setup();

  // get number of physical variables
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get maximum number of solution points in a cell
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  const CFuint nbrElemTypes = sdLocalData.size();
  CFuint maxNbrSolPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrSolPnts = sdLocalData[iElemType]->getNbrOfSolPnts();
    maxNbrSolPnts = maxNbrSolPnts > nbrSolPnts ? maxNbrSolPnts : nbrSolPnts;
  }

  // resize m_solPntJacobDets
  m_solPntJacobDets.resize(maxNbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::unsetup()
{
  CFAUTOTRACE;
  SpectralFDMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::execute()
{
  CFAUTOTRACE;

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  CFLogDebugMin("StdSourceTerm::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  // loop over elements
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = (*elemType)[m_iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[m_iElemType].getStartIdx();

    // get data needed for source term computation
    getSourceTermData();

    // loop over elements
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      // build the GeometricEntity
      geoData.idx = cellIdx;
      m_cell = geoBuilder->buildGE();

      // get the states
      m_cellStates = m_cell->getStates();

      // get jacobian determinants at solution points
      m_solPntJacobDets =
          m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

      // add the source term if the current cell is parallel updatable
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        addSourceTerm();
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSourceTerm::getSourceTermData()
{
    // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

    // get solution point local coordinates
  m_solPntsLocalCoords = sdLocalData[m_iElemType]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    StdSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
