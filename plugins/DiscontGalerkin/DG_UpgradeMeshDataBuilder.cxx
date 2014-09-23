#include "Common/PE.hh"

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/MeshData.hh"

#include "Framework/PhysicalModel.hh"

#include "DiscontGalerkin/DiscontGalerkin.hh"
#include "DiscontGalerkin/DG_UpgradeMeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////


Environment::ObjectProvider<DG_UpgradeMeshDataBuilder,
       MeshDataBuilder,
       DiscontGalerkinModule,
       1>
dG_meshUpgradeBuilderProvider("DG_MeshUpgrade");

//////////////////////////////////////////////////////////////////////////////

void DG_UpgradeMeshDataBuilder::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< std::string >("SolutionPolyOrder","Polynomial order for shape function in the discontinuous Galerkin aproximation.");
}

//////////////////////////////////////////////////////////////////////////////

DG_UpgradeMeshDataBuilder::DG_UpgradeMeshDataBuilder(const std::string& name) :
 DG_MeshDataBuilder(name),
 m_PolyOrder()
{
 addConfigOptionsTo(this);

 m_PolyOrderStr = "P2";
 setParameter( "SolutionPolyOrder", &m_PolyOrderStr );
}

//////////////////////////////////////////////////////////////////////////////

DG_UpgradeMeshDataBuilder::~DG_UpgradeMeshDataBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void DG_UpgradeMeshDataBuilder::releaseMemory()
{
 DG_MeshDataBuilder::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void DG_UpgradeMeshDataBuilder::configure ( Config::ConfigArgs& args )
{

 DG_MeshDataBuilder::configure(args);

 m_PolyOrder = CFPolyOrder::Convert::to_enum( m_PolyOrderStr );

}

//////////////////////////////////////////////////////////////////////////////

CFuint DG_UpgradeMeshDataBuilder::getNbrOfStatesInType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder)
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

void DG_UpgradeMeshDataBuilder::computeGeoTypeInfo()
{
 CFAUTOTRACE;

 // get the element type data
 SafePtr< std::vector< ElementTypeData > > elementType = getCFmeshData().getElementTypeData();

 std::vector< ElementTypeData >::iterator type_itr = elementType->begin();
 for (; type_itr != elementType->end(); ++type_itr)
 {
  // get the number of control volumes (states) in this element type
  const CFuint nbStatesPerElem = getNbrOfStatesInType(type_itr->getGeoShape(),m_PolyOrder);
  type_itr->setNbStates(nbStatesPerElem);
  type_itr->setSolOrder(m_PolyOrder);
 }

 /// @todo since only meshes with the same order are currently supported
 ///   we have to change the order in the CFmeshData and
 ///   not only on the ElementTypeData
 getCFmeshData().setSolutionPolyType(CFPolyForm::LAGRANGE);
 getCFmeshData().setSolutionPolyOrder(m_PolyOrder);

 // continue with the standard algorithm
 DG_MeshDataBuilder::computeGeoTypeInfo();
}

//////////////////////////////////////////////////////////////////////////////

void DG_UpgradeMeshDataBuilder::createTopologicalRegionSets()
{
 CFAUTOTRACE;

 // first transform the cell-states connectivity
 // from cell centered to a DiscontGalerkin
 upgradeStateConnectivity();

 // recreate the states
 recreateStates();

 // continue with the standart algorithm
 DG_MeshDataBuilder::createTopologicalRegionSets();
}

//////////////////////////////////////////////////////////////////////////////

void DG_UpgradeMeshDataBuilder::upgradeStateConnectivity()
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
 const CFuint nbStatesPerElem = getNbrOfStatesInType(type_itr->getGeoShape(),m_PolyOrder);

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
 const CFuint nbStatesPerElem = getNbrOfStatesInType(type_itr->getGeoShape(),m_PolyOrder);

  // get the number of elements
 const CFuint nbrElems = type_itr->getNbElems();

  // get start index of this element type in global element list
 CFuint globalIdx = type_itr->getStartIdx();

  // loop over elements of this type
 for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++globalIdx)
  {
   for (CFuint jState = 0; jState < nbStatesPerElem; ++jState, ++stateID)
   {
// CFout << "\n*********" << (*cellStates)(globalIdx,jState) << "  " << CFendl;
    (*cellStates)(globalIdx,jState) = stateID;
// CFout << (*cellStates)(globalIdx,jState) << CFendl;
   }
  }
 }

 cf_assert(stateID == columnPattern.sum());

 MeshDataStack::getActive()->setTotalStateCount(stateID);
}

//////////////////////////////////////////////////////////////////////////////

void DG_UpgradeMeshDataBuilder::recreateStates()
{
 CFAUTOTRACE;

 Common::SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

 DataHandle<State*,GLOBAL> states = getCFmeshData().getStatesHandle();

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
  states[iState]->setLocalID(iState);
  cf_assert(!Common::PE::GetPE().IsParallel());
  CFuint globalID = states[iState]->getLocalID();
  states[iState]->setGlobalID(globalID);
 }
}

//////////////////////////////////////////////////////////////////////////////

}// namespace DiscontGalerkin

}// namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
