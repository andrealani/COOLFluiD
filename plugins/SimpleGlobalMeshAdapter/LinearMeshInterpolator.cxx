// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/progress.hpp>

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "LinearMeshInterpolator.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/State.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntityPool.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider < LinearMeshInterpolator,
                        SimpleMeshAdapterData,
                        SimpleGlobalMeshAdapterModule>
aLinearMeshInterpolatorProvider("LinearMeshInterpolator");

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >("MinCoord","???");
  options.addConfigOption< std::vector<CFreal> >("MaxCoord","???");
  options.addConfigOption< std::vector<CFuint> >("NbSubdiv","Number of subdivisions");
}

//////////////////////////////////////////////////////////////////////////////

LinearMeshInterpolator::LinearMeshInterpolator(const std::string& name)  :
  SimpleMeshAdapterCom(name),
  socket_states("states"),
  socket_otherStates("states")
{

   addConfigOptionsTo(this);

   _minDomainCoord = std::vector<CFreal>();
   setParameter("MinCoord",&_minDomainCoord);

   _maxDomainCoord = std::vector<CFreal>();
   setParameter("MaxCoord",&_maxDomainCoord);

   _nbSubdiv = std::vector<CFuint>();
   setParameter("NbSubdiv",&_nbSubdiv);
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::configure ( Config::ConfigArgs& args )
{
  SimpleMeshAdapterCom::configure(args);

  socket_otherStates.setDataSocketNamespace(getMethodData().getOtherNamespace());
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::setup()
{
  CFAUTOTRACE;

  SimpleMeshAdapterCom::setup();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  cf_assert(_nbSubdiv.size() == nbDim);
  cf_assert(_minDomainCoord.size() == nbDim);
  cf_assert(_maxDomainCoord.size() == nbDim);

  for(CFuint iDim=0;iDim<nbDim;++iDim)
  {
    _minDomainCoord[iDim] -= MathTools::MathConsts::CFrealEps();
    _maxDomainCoord[iDim] += MathTools::MathConsts::CFrealEps();
  }
  _coord.resize(nbDim);
  _tempVector.resize(nbDim);
  _minCoord.resize(nbDim);
  _maxCoord.resize(nbDim);

  CFuint nbBoxes = _nbSubdiv[0];
  for(CFuint iDim = 1; iDim < nbDim ; ++iDim)
  {
    nbBoxes *= _nbSubdiv[iDim];
  }
  _boxes.resize(nbBoxes);
  _boxID.resize(nbDim);

  SetupBoxLimits();

  _interpolatedState.resize(nbEqs);

  std::string name = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance().getNamespace(name);
  const CFuint maxNbNodesInCell =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->Statistics().getMaxNbNodesInCell();

  _shapeFunctions.resize(maxNbNodesInCell);
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::execute()
{
  CFAUTOTRACE;

  CFout << "Interpolating solution on new mesh\n";
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();

  BoxClassification();

  boost::progress_display progressBar(nbStates);

  CFout << "State interpolation (" << nbStates << " states to interpolate)\n";
  for(CFuint iState=0; iState< nbStates; iState++)
  {

    ++(progressBar);
    _stateID = iState;
    _coord = states[iState]->getCoordinates();

    SearchInBoxes();
  }

}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::SetupBoxLimits()
{
  CFAUTOTRACE;

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  _boxesMin.resize(nbDim);
  _boxesMax.resize(nbDim);
  for(CFuint iDim = 0; iDim < nbDim ; ++iDim)
  {
    cf_assert(_nbSubdiv[iDim] > 0);
    cf_assert(_maxDomainCoord[iDim] > _minDomainCoord[iDim]);
    const CFreal interval = (_maxDomainCoord[iDim] - _minDomainCoord[iDim])/_nbSubdiv[iDim];

    _boxesMin[iDim].resize(_nbSubdiv[iDim]);
    _boxesMax[iDim].resize(_nbSubdiv[iDim]);
    for(CFuint iBox = 0; iBox < _nbSubdiv[iDim] ; ++iBox)
    {
      _boxesMin[iDim][iBox] = _minDomainCoord[iDim] + iBox*interval;
      _boxesMax[iDim][iBox] = _minDomainCoord[iDim] + (iBox+1)*interval;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::BoxClassification()
{
  CFAUTOTRACE;

  //Create a structured mesh dividing the domain into rectangular zones
  //which contains only a subset of cells
  //You do this once and it is then faster to search

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  std::string otherNamespace = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance().getNamespace(otherNamespace);

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->getTrs("InnerCells");

  const CFuint nbCells = cells->getLocalNbGeoEnts();
  GeometricEntityPool<StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setupInNamespace(otherNamespace);

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.trs = cells;

  const CFuint maxNbNodesInCell =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->Statistics().getMaxNbNodesInCell();
  std::vector<CFuint> cellBoxes(maxNbNodesInCell);

  //Loop over the cells in the whole domain and put the cells in the boxes
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    geoData.idx = iCell;
    GeometricEntity& cell = *geoBuilder.buildGE();
    const CFuint nbNodesInCell = cell.nbNodes();

//  CFout << "Cell: " << iCell <<"\n";

    //initialize to extreme values
    for(CFuint iDim = 0; iDim < nbDim; iDim++)
    {
      _minCoord[iDim] = MathTools::MathConsts::CFrealMax();
      _maxCoord[iDim] = -MathTools::MathConsts::CFrealMax();
    }

    //Compute the max and min coordinates
    for(CFuint iNode=0; iNode < nbNodesInCell; iNode++)
    {
      _coord = *(cell.getNode(iNode));
// CFout << "Coord: " << _coord <<"\n";
      for(CFuint iDim = 0; iDim < nbDim; iDim++)
      {
        _minCoord[iDim] = std::min(_minCoord[iDim], _coord[iDim]);
        _maxCoord[iDim] = std::max(_maxCoord[iDim], _coord[iDim]);
      }
    }

// CFout << "Min Coord: " << _minCoord <<"\n";
// CFout << "Max Coord: " << _maxCoord <<"\n";
// CFout << "Trying to find the boxes that contain the cell\n";

    findBoxesFromMinMax();

    for(CFuint iBox=0;iBox<_listBoxes.size(); ++iBox)
    {
      _boxes[_listBoxes[iBox]].push_back(iCell);
//  CFout << "Put in Box: " << _listBoxes[iBox] << "\n";
    }

    //release the GeometricEntity
    geoBuilder.releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::SearchInBoxes()
{

  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

// CFout << "Searching for state at coord: " << _coord <<"\n";
  const CFuint boxIndex = findBoxFromCoord();
  const CFuint nbCellsInBox = _boxes[boxIndex].size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  std::string name = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance().getNamespace(name);

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->getTrs("InnerCells");

  GeometricEntityPool<StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setupInNamespace(name);

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.trs = cells;

  //Loop over the cells in the box
  bool elementFound(false);
  for (CFuint iCell = 0; (iCell < nbCellsInBox) && (!elementFound); ++iCell) {

//   CFout << "NbCellsInBox: " << nbCellsInBox <<"\n";
    // build the GeometricEntity
    const CFuint cellID = _boxes[boxIndex][iCell];
    geoData.idx = cellID;
    GeometricEntity *const currCell = geoBuilder.buildGE();

    //check if coord is in cell
    elementFound = currCell->isInElement(_coord);
    if(elementFound)
    {
      _shapeFunctions = currCell->computeShapeFunctionAtCoord(_coord);

      std::vector<State*>* cellStates = currCell->getStates();
      const CFuint nbNodes = currCell->nbNodes();
      const CFuint nbStates = currCell->nbStates();

      //this is for FEM/RDS only!!!
      cf_assert(nbNodes == nbStates);
      for (CFuint j = 0; j < nbEqs; ++j)
      {
        _interpolatedState[j] = _shapeFunctions[0] * ((*(*cellStates)[0])[j]);
        for (CFuint k = 1; k < nbNodes; ++k)
        {
          _interpolatedState[j] += _shapeFunctions[k] * ((*((*cellStates)[k]))[j]);
        }
      }
      (states[_stateID])->copyData(_interpolatedState);
    }

    //release the GeometricEntity
    geoBuilder.releaseGE();
  }

  //if false: compute the closest point (projection)
  if(!elementFound){
    CFout << "No element was found for state: " <<_coord <<" in Box: "<< boxIndex<<"\n";

    for (CFuint iCell = 0; (iCell < nbCellsInBox); ++iCell) {
      CFout << "Containing cell: " << _boxes[boxIndex][iCell] <<"\n";
    }
    computeClosestValue();

    (states[_stateID])->copyData(_interpolatedState);
  }

}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::computeClosestValue()
{
//on which boundary should we look??
//Look on the boundaries with the same
cf_assert(false);

 DataHandle<State*,GLOBAL> otherStates = socket_otherStates.getDataHandle();

///@todo compute using the 2 (DIM_2D) or 3(DIM_3D) nearest states
///@todo compute using the projection on the nearest face

  std::string name = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance().getNamespace(name);

  std::vector< Common::SafePtr<TopologicalRegionSet> > trsList =
    MeshDataStack::getInstance().getEntryByNamespace(otherNsp)->getTrsList();

  const CFuint nbTrsInOtherMesh = trsList.size();

  CFreal minDistance = MathTools::MathConsts::CFrealMax();
  for (CFuint iTrs = 0; iTrs < nbTrsInOtherMesh; ++iTrs) {

    Common::SafePtr<TopologicalRegionSet> trs = trsList[iTrs];

    const std::string trsName = trs->getName();
    if((trsName != "InnerCells") && (trsName != "InnerCells")){

      const CFuint nbOtherStates = trs->getStatesInTrs()->size();
      for(CFuint iOtherState=0; iOtherState < nbOtherStates; iOtherState++)
      {
        const CFuint otherStateID = (*(trs->getStatesInTrs()))[iOtherState];
        _tempVector = (otherStates[otherStateID])->getCoordinates() - _coord;
        const CFreal distance = _tempVector.norm2();

        if(distance < minDistance)
        {
          _interpolatedState = *(otherStates[otherStateID]);
          minDistance = distance;
        }
      }
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

CFuint LinearMeshInterpolator::findBoxIndex()
{
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  //find the index using _boxID
  if(nbDim == DIM_2D)
  {
    return _boxID[0] + _boxID[1]*_nbSubdiv[0];
  }
  if(nbDim == DIM_3D)
  {
    return _boxID[0] + _boxID[1]*_nbSubdiv[0] + (_boxID[2]*_nbSubdiv[1]);
  }

  cf_assert(false);
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

CFuint LinearMeshInterpolator::findBoxFromCoord()
{
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

///@todo modify this to be faster...
/// should be possible to compute _boxID[iDim] for each dimension based on:
/// _nbSubdiv[iDim]
/// _minDomainCoord[iDim]
/// _maxDomainCoord[iDim]

  for(CFuint iDim = 0; iDim < nbDim ; ++iDim)
  {
    cf_assert(_nbSubdiv[iDim] > 0);
    const CFreal interval = (_maxDomainCoord[iDim] - _minDomainCoord[iDim])/_nbSubdiv[iDim];

    _boxID[iDim] = static_cast<CFuint>(std::floor((_coord[iDim] - _minDomainCoord[iDim])/interval));
  }
/*
  for(CFuint iDim = 0; iDim < nbDim ; ++iDim)
  {
    bool boxFound = false;
    for(CFuint iBox = 0; iBox < _nbSubdiv[iDim] ; ++iBox)
    {
      const CFreal minBox = _boxesMin[iDim][iBox];
      const CFreal maxBox = _boxesMax[iDim][iBox];

      if((_coord[iDim] >= minBox) && (_coord[iDim] <= maxBox))
      {
        _boxID[iDim] = iBox;
        boxFound = true;
      }
    }
    cf_assert(boxFound);
  }*/

  return findBoxIndex();
}

//////////////////////////////////////////////////////////////////////////////

void LinearMeshInterpolator::findBoxesFromMinMax()
{
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  vector<CFuint> minBoxID(nbDim);
  vector<CFuint> maxBoxID(nbDim);

  CFuint totalNbBoxes = 1;
  for(CFuint iDim = 0; iDim < nbDim ; ++iDim)
  {
    cf_assert(_nbSubdiv[iDim] > 0);
    cf_assert(_maxCoord[iDim] > _minCoord[iDim]);
    const CFreal interval = (_maxDomainCoord[iDim] - _minDomainCoord[iDim])/_nbSubdiv[iDim];

    minBoxID[iDim] = static_cast<CFuint>(std::floor((_minCoord[iDim] - _minDomainCoord[iDim])/interval));
    maxBoxID[iDim] = static_cast<CFuint>(std::ceil((_maxCoord[iDim] - _minDomainCoord[iDim])/interval));

    totalNbBoxes *= (maxBoxID[iDim] - minBoxID[iDim]);
  }

  //Compute the array of boxes to use
  _listBoxes.resize(totalNbBoxes);

  CFuint index = 0;
  for(CFuint iX = minBoxID[0]; iX < maxBoxID[0] ; ++iX)
  {
    _boxID[0] = iX;
    for(CFuint iY = minBoxID[1]; iY < maxBoxID[1] ; ++iY)
    {
      _boxID[1] = iY;
      if(nbDim == DIM_2D)
      {
        _listBoxes[index] = findBoxIndex();
        ++index;
      }
      if(nbDim == DIM_3D)
      {
        for(CFuint iZ = minBoxID[2]; iZ < maxBoxID[2] ; ++iZ)
        {
          _boxID[2] = iZ;
          _listBoxes[index] = findBoxIndex();
          ++index;
        }
      }
    }
  }


}


//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  LinearMeshInterpolator::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_states);
  result.push_back(&socket_otherStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
