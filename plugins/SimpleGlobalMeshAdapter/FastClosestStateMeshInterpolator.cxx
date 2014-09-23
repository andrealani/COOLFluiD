// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "MathTools/MathConsts.hh"
#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"
#include "SimpleGlobalMeshAdapter/FastClosestStateMeshInterpolator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FastClosestStateMeshInterpolator,
		      SimpleMeshAdapterData,
		      SimpleGlobalMeshAdapterModule>
FastClosestStateMeshInterpolatorProvider("FastClosestStateMeshInterpolator");

//////////////////////////////////////////////////////////////////////////////

void FastClosestStateMeshInterpolator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >("MinCoord","Name of the journal file.");
  options.addConfigOption< std::vector<CFreal> >("MaxCoord","Name of the output mesh.");
  options.addConfigOption< std::vector<CFuint> >("NbSubdiv","Name of the output mesh.");
}

//////////////////////////////////////////////////////////////////////////////

FastClosestStateMeshInterpolator::FastClosestStateMeshInterpolator(const std::string& name)  :
  SimpleMeshAdapterCom(name),
  socket_states("states"),
  socket_otherStates("states"),
  _dimension(0)
{
  addConfigOptionsTo(this);

  _minCoord = std::vector<CFreal>();
  setParameter("MinCoord",&_minCoord);

  _maxCoord = std::vector<CFreal>();
  setParameter("MaxCoord",&_maxCoord);

  _nbSubdiv = std::vector<CFuint>();
  setParameter("NbSubdiv",&_nbSubdiv);
}

//////////////////////////////////////////////////////////////////////////////

void FastClosestStateMeshInterpolator::configure ( Config::ConfigArgs& args )
{
  SimpleMeshAdapterCom::configure(args);

  socket_otherStates.setDataSocketNamespace(getMethodData().getOtherNamespace());
}

//////////////////////////////////////////////////////////////////////////////

void FastClosestStateMeshInterpolator::setup()
{
  CFAUTOTRACE;

  SimpleMeshAdapterCom::setup();

  _dimension = PhysicalModelStack::getActive()->getDim();
  cf_assert(_dimension > DIM_1D);
  cf_assert(_nbSubdiv.size() == _dimension);
  cf_assert(_minCoord.size() == _dimension);
  cf_assert(_maxCoord.size() == _dimension);

  _coord.resize(_dimension );
  _tempVector.resize(_dimension );

  CFuint nbBoxes = _nbSubdiv[0];
  for(CFuint iDim = 1; iDim < _dimension  ; ++iDim)
  {
    nbBoxes *= _nbSubdiv[iDim];
  }
  _boxes.resize(nbBoxes);
  _boxID.resize(_dimension);
  _boxIDTol.resize(_dimension);

  SetupBoxLimits();
}

//////////////////////////////////////////////////////////////////////////////

void FastClosestStateMeshInterpolator::execute()
{
  CFAUTOTRACE;

  CFout << "Interpolating solution on new mesh\n";

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();

  CFout << "Preselection in boxes\n";
  BoxClassification();

  CFout << "States interpolation\n";
  for(CFuint iState=0; iState< nbStates; iState++)
  {
    //  CFout << iState << "\n";
    _stateID = iState;
    _coord = states[iState]->getCoordinates();

    SearchInBoxes();
  }

}

//////////////////////////////////////////////////////////////////////////////

void FastClosestStateMeshInterpolator::SetupBoxLimits()
{
  _boxesMin.resize(_dimension);
  _boxesMax.resize(_dimension);
  for(CFuint iDim = 0; iDim < _dimension ; ++iDim)
  {
    cf_assert(_nbSubdiv[iDim] > 0);
    cf_assert(_maxCoord[iDim] > _minCoord[iDim]);
    const CFreal interval = (_maxCoord[iDim] - _minCoord[iDim])/_nbSubdiv[iDim];

    _boxesMin[iDim].resize(_nbSubdiv[iDim]);
    _boxesMax[iDim].resize(_nbSubdiv[iDim]);
    for(CFuint iBox = 0; iBox < _nbSubdiv[iDim] ; ++iBox)
    {
      _boxesMin[iDim][iBox] = _minCoord[iDim] + iBox*interval;
      _boxesMax[iDim][iBox] = _minCoord[iDim] + (iBox+1)*interval;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FastClosestStateMeshInterpolator::BoxClassification()
{
  //Create a structured mesh dividing the domain into rectangular zones

  DataHandle<State*,GLOBAL> otherStates = socket_otherStates.getDataHandle();
  const CFuint nbOtherStates = otherStates.size();

  for(CFuint iOtherState=0; iOtherState < nbOtherStates; iOtherState++)
  {
    _coord = (otherStates[iOtherState])->getCoordinates();
    findBoxFromCoordTolerance();

    for(CFuint i=0; i < _belongToBoxesID.size(); i++)
    {
      (_boxes[_belongToBoxesID[i]]).push_back(iOtherState);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FastClosestStateMeshInterpolator::SearchInBoxes()
{
  DataHandle<State*,GLOBAL> otherStates = socket_otherStates.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint boxIndex = findBoxFromCoord();
  const CFuint nbStatesInBox = _boxes[boxIndex].size();

  //Loop over the states in the box
  CFreal minDistance = MathTools::MathConsts::CFrealMax();
  for(CFuint iState=0; iState < nbStatesInBox; ++iState)
  {
    const CFuint otherStateID = _boxes[boxIndex][iState];
    _tempVector = _coord - (otherStates[otherStateID])->getCoordinates();
    const CFreal distance = _tempVector.norm2();

    if(distance < minDistance){
      minDistance = distance;
      (states[_stateID])->copyData(*(otherStates[otherStateID]));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

///@this is to allow a state to belong to several boxes...
void FastClosestStateMeshInterpolator::findBoxFromCoordTolerance()
{
  _belongToBoxesID.resize(0);

  for(CFuint iDim = 0; iDim < _dimension ; ++iDim)
  {
    const CFuint nbSubDim = _nbSubdiv[iDim];
    const CFreal interval = (_maxCoord[iDim] - _minCoord[iDim])/nbSubDim;
    const CFreal tolerance = 0.5 * interval;

    _boxIDTol[iDim].resize(0);
    for(CFuint iBox = 0; iBox < nbSubDim; ++iBox)
    {
      const CFreal minBox = _boxesMin[iDim][iBox] - tolerance;
      const CFreal maxBox = _boxesMax[iDim][iBox] + tolerance;

      if((_coord[iDim] > minBox) && (_coord[iDim] <= maxBox))
      {
        _boxIDTol[iDim].push_back(iBox);
      }
    }
  }

  if(_dimension == DIM_3D)
  {
    for(CFuint iBox = 0; iBox < _boxIDTol[0].size(); ++iBox)
    {
      for(CFuint jBox = 0; jBox < _boxIDTol[1].size(); ++jBox)
      {
        for(CFuint kBox = 0; kBox < _boxIDTol[2].size(); ++kBox)
        {
          _boxID[0] = _boxIDTol[0][iBox];
          _boxID[1] = _boxIDTol[1][jBox];
          _boxID[2] = _boxIDTol[2][kBox];
          _belongToBoxesID.push_back(findBoxIndex());
        }
      }
    }
  }
  else {
    cf_assert(_dimension == DIM_2D);

    for(CFuint iBox = 0; iBox < _boxIDTol[0].size(); ++iBox)
    {
      for(CFuint jBox = 0; jBox < _boxIDTol[1].size(); ++jBox)
      {
	_boxID[0] = _boxIDTol[0][iBox];
	_boxID[1] = _boxIDTol[1][jBox];
	_belongToBoxesID.push_back(findBoxIndex());
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  FastClosestStateMeshInterpolator::needsSockets()
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
