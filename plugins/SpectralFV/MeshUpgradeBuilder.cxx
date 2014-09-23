#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/MeshUpgradeBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////


Environment::ObjectProvider<MeshUpgradeBuilder,
               MeshDataBuilder,
               SpectralFVModule,
               1>
meshUpgradeBuilderProvider("MeshUpgrade");

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SVPolynomialOrder","Spectral finite volume polynomial order.");
}

//////////////////////////////////////////////////////////////////////////////

MeshUpgradeBuilder::MeshUpgradeBuilder(const std::string& name) :
  SpectralFVBuilder(name),
  m_svPolyOrder()
{
  addConfigOptionsTo(this);

  m_svPolyOrderStr = "P1";
  setParameter( "SVPolynomialOrder", &m_svPolyOrderStr );
}

//////////////////////////////////////////////////////////////////////////////

MeshUpgradeBuilder::~MeshUpgradeBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::releaseMemory()
{
  SpectralFVBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::configure ( Config::ConfigArgs& args )
{
  SpectralFVBuilder::configure(args);

  m_svPolyOrder = CFPolyOrder::Convert::to_enum( m_svPolyOrderStr );
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshUpgradeBuilder::getNbrOfStatesInSVType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder)
{
  switch (geoShape)
  {
    case CFGeoShape::LINE:
      return polyOrder + 1;
    case CFGeoShape::TRIAG:
      return (polyOrder + 1)*(polyOrder + 2)/2;
    case CFGeoShape::QUAD:
      return (polyOrder + 1)*(polyOrder + 1);
    case CFGeoShape::TETRA:
      return (polyOrder + 1)*(polyOrder + 2)*(polyOrder + 3)/6;
    case CFGeoShape::HEXA:
      return (polyOrder + 1)*(polyOrder + 1)*(polyOrder + 1);
    default:
      throw Common::NotImplementedException (FromHere(),"Element type not implemented...");
  }
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::computeGeoTypeInfo()
{
  CFAUTOTRACE;

  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();

  std::vector< ElementTypeData >::iterator type_itr = elementType->begin();
  for (; type_itr != elementType->end(); ++type_itr)
  {
    // get the number of control volumes (states) in this element type
    const CFuint nbStatesPerElem = getNbrOfStatesInSVType(type_itr->getGeoShape(),m_svPolyOrder);
    type_itr->setNbStates(nbStatesPerElem);
    type_itr->setSolOrder(m_svPolyOrder);
  }

  /// @todo since only meshes with the same order are currently supported
  ///       we have to change the order in the CFmeshData and
  ///       not only on the ElementTypeData
  getCFmeshData().setSolutionPolyOrder(m_svPolyOrder);

  // continue with the standard algorithm
  SpectralFVBuilder::computeGeoTypeInfo();
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::createTopologicalRegionSets()
{
  CFAUTOTRACE;

  // first transform the cell-states connectivity
  // from cell centered to a SpectralFV
  upgradeStateConnectivity();

  // recreate the states
  recreateStates();

  // continue with the standart algorithm
  SpectralFVBuilder::createTopologicalRegionSets();
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::upgradeStateConnectivity()
{
  CFAUTOTRACE;

  // get the element type data
  SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();
  std::vector< ElementTypeData >::iterator type_itr;

  // get the connectivity that we will change
  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbElems = getNbElements();

  // we will create a mesh with equal order on all the elements
  // and we reset the connectivity accordingly, using a std::valarray
  // first create the std::valarray
  std::valarray<CFuint> columnPattern(nbElems);
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of control volumes (states) in this element type
    const CFuint nbStatesPerElem = getNbrOfStatesInSVType(type_itr->getGeoShape(),m_svPolyOrder);

    // get the number of elements
    const CFuint nbrElems = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
    {
      columnPattern[globalIdx] = nbStatesPerElem;
    }
  }
  // then resize the connectivity
  cellStates->resize(columnPattern);

  // loop on all the elements and reassign the IDs of the states
  CFuint stateID = 0;
  for (type_itr = elementType->begin(); type_itr != elementType->end(); ++type_itr)
  {
    // get the number of control volumes (states) in this element type
    const CFuint nbStatesPerElem = getNbrOfStatesInSVType(type_itr->getGeoShape(),m_svPolyOrder);

    // get the number of elements
    const CFuint nbrElems = type_itr->getNbElems();

    // get start index of this element type in global element list
    CFuint globalIdx = type_itr->getStartIdx();

    // loop over elements of this type
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
    {
      for (CFuint jState = 0; jState < nbStatesPerElem; ++jState, ++stateID)
      {
        (*cellStates)(globalIdx,jState) = stateID;
      }
    }
  }

  cf_assert(stateID == columnPattern.sum());
}

//////////////////////////////////////////////////////////////////////////////

void MeshUpgradeBuilder::recreateStates()
{
  CFAUTOTRACE;

  Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  DataHandle < Framework::State*, Framework::GLOBAL > states = getCFmeshData().getStatesHandle();

  const CFuint newNbStates = cellStates->size();
  const CFuint oldNbStates = states.size();

  // delete the existing states
  for (CFuint i = 0; i < oldNbStates; ++i)
  {
    deletePtr(states[i]);
  }
  IndexList<State>::getList().reset();

  // Resize the datahandle for the states
  states.resize(newNbStates);

  // allocate the new states
  const CFuint nbeq = PhysicalModelStack::getActive()->getNbEq();
  RealVector stateData (nbeq);
  for (CFuint iState = 0; iState < states.size(); ++iState)
  {
    getCFmeshData().createState(iState,stateData);
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

