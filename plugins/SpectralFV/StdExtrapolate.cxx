#include "SpectralFV/StdExtrapolate.hh"
#include "SpectralFV/SpectralFV.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdExtrapolate,SpectralFVMethodData,SpectralFVModule >
  StdExtrapolateProvider("StdExtrapolate");

//////////////////////////////////////////////////////////////////////////////

StdExtrapolate::StdExtrapolate(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_nstates("nstates"),
  socket_nstatesProxy("nstatesProxy"),
  socket_states("states"),
  socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdExtrapolate::~StdExtrapolate()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StdExtrapolate::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_nstatesProxy);
  result.push_back(&socket_nstates);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  StdExtrapolate::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdExtrapolate::execute()
{
  CFAUTOTRACE;

  // get the datahandle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get the datahandle for the nodal states
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();

  // get cell-node connectivity
  SafePtr<MeshData::ConnTable> cellNodes = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");

  // get cell-state connectivity
  SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the number of nodes
  const CFuint nbrNodes = nstates.size();

  // variable to count the number of cells that contribute to the state in a node
  /// @todo should this variable be a class variable, to avoid repeated memory allocation and release
  /// during run? How often is this function called?
  vector< CFuint > nbrNodeCells(nbrNodes,0);

  // reset the nodal states to zero
  for (CFuint iNState = 0; iNState < nbrNodes; ++iNState)
  {
    nstates[iNState] = 0.0;
  }

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  // only element types with a P1 geometry at most can be treated in this loop
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get SV nodes local coordinates
    SafePtr< vector< RealVector > > svNodeLocalCoords = svLocalData[iElemType]->getSVNodeCoords();

    // number of nodes to this element type
    const CFuint nbrCellNodes = svNodeLocalCoords->size();

    // get SV polynomial values at SV nodes
    vector< vector< CFreal > > svPolyValsAtSVNodes
                                  = svLocalData[iElemType]->getSVPolyValsAtNode(*svNodeLocalCoords);

    // get the number of control volumes per spectral volume
    const CFuint nbrCVs = svLocalData[iElemType]->getNbrOfCVs();

    // get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      const CFuint stateStartIdx = (*cellStates)(cellIdx,0); // states per element are ordered sequentially

      for (CFuint iNode = 0; iNode < nbrCellNodes; ++iNode)
      {
        // get node local index
        const CFuint nodeIdx = (*cellNodes)(cellIdx,iNode);

        // increase nbrNodeCells
        ++nbrNodeCells[nodeIdx];

        // add contribution of cell to this node
        for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
        {
          nstates[nodeIdx] += svPolyValsAtSVNodes[iNode][iCV]*(*states[stateStartIdx+iCV]);
        }
      }
    }
  }

  // divide each nodal state by the number of cells that contributed to it
  for (CFuint iNState = 0; iNState < nbrNodes; ++iNState)
  {
    nstates[iNState] /= nbrNodeCells[iNState];
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD
