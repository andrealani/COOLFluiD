// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <numeric>
#include <boost/progress.hpp>

#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"

#include "Common/NoSuchValueException.hh"

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"
#include "SimpleGlobalMeshAdapter/ShepardMeshInterpolator.hh"
#include "MathTools/MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ShepardMeshInterpolator,
		      SimpleMeshAdapterData,
		      SimpleGlobalMeshAdapterModule>
ShepardMeshInterpolatorProvider("ShepardMeshInterpolator");

//////////////////////////////////////////////////////////////////////////////

void ShepardMeshInterpolator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFuint> >("NbSubdiv","Number of subdivision per dimension.");
  options.addConfigOption< CFuint >("NbSelectedStates","Number of states which contribute to the interpolation.");
  options.addConfigOption< std::vector<CFreal> >("MinCoord","Name of the journal file.");
  options.addConfigOption< std::vector<CFreal> >("MaxCoord","Name of the output mesh.");
}

//////////////////////////////////////////////////////////////////////////////

ShepardMeshInterpolator::ShepardMeshInterpolator(const std::string& name)  :
  SimpleMeshAdapterCom(name),
  socket_states("states"),
  socket_otherStates("states"),
  _tempVector(),
  _minMaxCoordInBox(),
  _stateIDsInBox(),
  _globalBoxMinMax()
{
  addConfigOptionsTo(this);

  _nbSubdiv = std::vector<CFuint>();
  setParameter("NbSubdiv",&_nbSubdiv);

  _nbSelectedStates = 5;
  setParameter("NbSelectedStates",&_nbSelectedStates);

  _minCoord = std::vector<CFreal>();
  setParameter("MinCoord",&_minCoord);

  _maxCoord = std::vector<CFreal>();
  setParameter("MaxCoord",&_maxCoord);
}

//////////////////////////////////////////////////////////////////////////////

void ShepardMeshInterpolator::configure ( Config::ConfigArgs& args )
{
  SimpleMeshAdapterCom::configure(args);

  socket_otherStates.setDataSocketNamespace(getMethodData().getOtherNamespace());
}

//////////////////////////////////////////////////////////////////////////////

void ShepardMeshInterpolator::setup()
{
  CFAUTOTRACE;

  SimpleMeshAdapterCom::setup();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(_nbSubdiv.size() == dim);
  _tempVector.resize(dim);

  setupBoxes();
}

//////////////////////////////////////////////////////////////////////////////

void ShepardMeshInterpolator::execute()
{
  CFAUTOTRACE;

  CFout << "Interpolating solution on new mesh\n";

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
//   const CFuint nbStates = states.size();

  CFout << "Preselection in boxes\n";
  BoxClassification();

  CFout << "States interpolation\n";
  SearchInBoxes();
}

//////////////////////////////////////////////////////////////////////////////

void ShepardMeshInterpolator::setupBoxes()
{
  // states in the old mesh
  DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle();

  // states in the new mesh
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  CFuint nbBoxes = _nbSubdiv[0];
  for(CFuint iDim = 1; iDim < dim ; ++iDim) {
    nbBoxes *= _nbSubdiv[iDim];
  }

  // starting guess for the number of old mesh states associated to one box
  // this allows to use vector::push_back more efficiently
  const CFuint avgNbStatesInBox = otherStates.size()/nbBoxes;
  _minMaxCoordInBox.resize(nbBoxes*dim*2);

  // storage associating each old stateID with a box ID
  _stateIDsInBox.resize(nbBoxes);
  for(CFuint i = 0; i < nbBoxes; ++i) {
    _stateIDsInBox[i].reserve(avgNbStatesInBox);
  }

  // set the limits of the bounding box containing the whole domain
  _globalBoxMinMax.resize(dim,2);
  for(CFuint iDim = 0; iDim < dim ; ++iDim) {
    _globalBoxMinMax(iDim,0) = _minCoord[iDim];
    _globalBoxMinMax(iDim,1) = _maxCoord[iDim];
  }

  // apply a 1% tolerance on the bounding box to be sure to include all nodes
  const CFreal factor = 1.01;
  RealVector interval(dim);
  
  for(CFuint iDim = 0; iDim < dim ; ++iDim) {
    (_globalBoxMinMax(iDim,0) <= 0.0) ?
      _globalBoxMinMax(iDim,0) *= factor : _globalBoxMinMax(iDim,0) /= factor;
    
    (_globalBoxMinMax(iDim,1) <= 0.0) ?
      _globalBoxMinMax(iDim,1) /= factor : _globalBoxMinMax(iDim,1) *= factor;
    
    cf_assert(_globalBoxMinMax(iDim,0) < _globalBoxMinMax(iDim,1));
    
    interval[iDim] = (_globalBoxMinMax(iDim,1) - _globalBoxMinMax(iDim,0))/_nbSubdiv[iDim];
  }
  
  CFLogNotice("ShepardMeshInterpolator::setupBoxes() => bounding box [min_Xi, max_Xi] = " << "\n" << _globalBoxMinMax << "\n");
  CFLogNotice("ShepardMeshInterpolator::setupBoxes() => interval [x y z] = " << interval << "\n");

  // ----------------- SANITY CHECKS ------------------------------ //
  // scan otherStates for checking if the bounding box limits are ok
  for (CFuint i = 0; i < otherStates.size(); ++i) {
    RealVector& node = otherStates[i]->getCoordinates();
    for(CFuint iDim = 0; iDim < dim ; ++iDim) {
      if (node[iDim] < _globalBoxMinMax(iDim,0)){
	throw NoSuchValueException(FromHere(),"node coordinates out of bounding box (min)");
      }
      if (node[iDim] > _globalBoxMinMax(iDim,1)){
	throw NoSuchValueException(FromHere(),"node coordinates out of bounding box (max)");
      }
    }
  }
  
  // scan states for checking if the bounding box limits are ok
  for (CFuint i = 0; i < states.size(); ++i) {
    RealVector& node = states[i]->getCoordinates();
    for(CFuint iDim = 0; iDim < dim ; ++iDim) {
      if (node[iDim] < _globalBoxMinMax(iDim,0)){
	throw NoSuchValueException(FromHere(),"node coordinates out of bounding box (min)");
      }
      if (node[iDim] > _globalBoxMinMax(iDim,1)){
	throw NoSuchValueException(FromHere(),"node coordinates out of bounding box (max)");
      }
    }
  }
  // -------------------------------------------------------------- //
  
  // define the boxes limits
  const CFuint dim2 = dim*2;
  RealMatrix boxLimits(dim, 2, &_minMaxCoordInBox[0]);
  
  if (dim == DIM_2D) {
    const CFuint nbSubX = _nbSubdiv[XX];
    const CFuint nbSubY = _nbSubdiv[YY];
    for (CFuint i = 0; i < nbSubX; ++i) {
      for (CFuint j = 0; j < nbSubY; ++j) {
	const CFuint globalBoxID = j*nbSubX + i;
	boxLimits.wrap(dim, 2, &_minMaxCoordInBox[globalBoxID*dim2]);
	boxLimits(XX,0) = _globalBoxMinMax(XX,0) + i*interval[XX];
	boxLimits(XX,1) = _globalBoxMinMax(XX,0) + (i+1)*interval[XX];
	boxLimits(YY,0) = _globalBoxMinMax(YY,0) + j*interval[YY];
	boxLimits(YY,1) = _globalBoxMinMax(YY,0) + (j+1)*interval[YY];
      }
    }
  }

  if (dim == DIM_3D) {
    const CFuint nbSubX = _nbSubdiv[XX];
    const CFuint nbSubY = _nbSubdiv[YY];
    const CFuint nbSubZ = _nbSubdiv[ZZ];

    for (CFuint i = 0; i < nbSubX; ++i) {
      for (CFuint j = 0; j < nbSubY; ++j) {
	for (CFuint k = 0; k < nbSubZ; ++k) {
	  const CFuint globalBoxID = k*nbSubX*nbSubY + j*nbSubX + i;
	  boxLimits.wrap(dim, 2, &_minMaxCoordInBox[globalBoxID*dim2]);
	  boxLimits(XX,0) = _globalBoxMinMax(XX,0) + i*interval[XX];
	  boxLimits(XX,1) = _globalBoxMinMax(XX,0) + (i+1)*interval[XX];
	  boxLimits(YY,0) = _globalBoxMinMax(YY,0) + j*interval[YY];
	  boxLimits(YY,1) = _globalBoxMinMax(YY,0) + (j+1)*interval[YY];
	  boxLimits(ZZ,0) = _globalBoxMinMax(ZZ,0) + k*interval[ZZ];
	  boxLimits(ZZ,1) = _globalBoxMinMax(ZZ,0) + (k+1)*interval[ZZ];
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ShepardMeshInterpolator::BoxClassification()
{
  // assign each one of the existing old states to a box
  DataHandle<State*,GLOBAL> otherStates = socket_otherStates.getDataHandle();
  
  const CFuint nbOtherStates = otherStates.size();
  const CFuint nbBoxes = _stateIDsInBox.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint dim2 = dim*2;
  RealMatrix boxLimits(dim, 2, &_minMaxCoordInBox[0]);
  
  CFuint countStateIDs = 0;
  for(CFuint iOtherState = 0; iOtherState < nbOtherStates; iOtherState++) {
    RealVector& coord = otherStates[iOtherState]->getCoordinates();

    bool foundBox = false;
    for(CFuint iBox = 0; iBox < nbBoxes; ++iBox) {
      boxLimits.wrap(dim, 2, &_minMaxCoordInBox[iBox*dim2]);
      CFuint countDimMatching = 0;
      for(CFuint iDim = 0; iDim < dim ; ++iDim) {
	if((coord[iDim] >=  boxLimits(iDim,0)) && (coord[iDim] <=  boxLimits(iDim,1))) {
	  countDimMatching++;
	}
      }

      // if the given state satisfies all coordinates constraints
      if (countDimMatching == dim) {
	_stateIDsInBox[iBox].push_back(iOtherState);
	countStateIDs++;
	foundBox = true;
	break;
      }
    }
    cf_assert(foundBox);
  }

  cf_always_assert(countStateIDs == nbOtherStates);
}

//////////////////////////////////////////////////////////////////////////////

void ShepardMeshInterpolator::SearchInBoxes()
{
  DataHandle<State*,GLOBAL> otherStates = socket_otherStates.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbBoxes = _stateIDsInBox.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint dim2 = dim*2;
  RealMatrix boxLimits(dim, 2, &_minMaxCoordInBox[0]);

  CFuint maxNbStatesInBox = 0;
  for(CFuint iBox = 0; iBox < nbBoxes; ++iBox) {
    maxNbStatesInBox = std::max(maxNbStatesInBox, static_cast<CFuint>(_stateIDsInBox[iBox].size()));
  }

  vector< pair<CFreal,State*> > distances;
  distances.reserve(maxNbStatesInBox);

  pair<CFreal, State*> tmpPair;
  CFuint countStateIDs = 0;
  
  boost::progress_display progressBar (nbStates);
  for(CFuint i = 0; i < nbStates; ++i) {
    
    ++(progressBar);
    RealVector& coord = states[i]->getCoordinates();

    //cout << "new coord[" << i << "] = " << coord << endl;
    //cout << "matching coord: ";
    bool foundBox = false;
    for(CFuint iBox = 0; iBox < nbBoxes; ++iBox) {
      boxLimits.wrap(dim, 2, &_minMaxCoordInBox[iBox*dim2]);
      CFuint countDimMatching = 0;

      for(CFuint iDim = 0; iDim < dim ; ++iDim) {
	if((coord[iDim] >=  boxLimits(iDim,0)) && (coord[iDim] <=  boxLimits(iDim,1))) {
	  countDimMatching++;
	}
      }

      // if the given state belongs to the current box (it satisfies all coordinates constraints)
      if (countDimMatching == dim) {
	const CFuint nbOtherStatesInBox = _stateIDsInBox[iBox].size();
	cf_assert(nbOtherStatesInBox  > 0);

	// this existing state has the exact same coordinates of the current one
	State* twinState = CFNULL;

	for (CFuint iState = 0; iState < nbOtherStatesInBox; ++iState) {
	  const CFuint otherStateID = _stateIDsInBox[iBox][iState];
	  State* const currOtherState = otherStates[otherStateID];
	  RealVector& otherCoord = currOtherState->getCoordinates();

	  //cout << "iBox = " << iBox << " => (" << otherCoord << "), ";
	  tmpPair.first  = MathFunctions::getDistance(coord, otherCoord);

	  if (!MathChecks::isZero(tmpPair.first)) {
	    tmpPair.second = currOtherState;
	    distances.push_back(tmpPair);
	  }
	  else {
	    twinState = currOtherState;
	  }
	}

	if (twinState == CFNULL) {
	  // sort the distances and take the closest _nbSelectedStates
	  std::sort(distances.begin(),distances.end());

	  CFreal sumWeights = 0.0;
	  *states[i] = 0.0;

	  //cout << "nbOtherStatesInBox = " << nbOtherStatesInBox << endl;
	  const CFuint selectionSize = std::min(nbOtherStatesInBox, _nbSelectedStates);

	  // take into account only the first selectionSize states != current one
	  CFuint counter = 0;
	  for (CFuint iState = 0; iState < nbOtherStatesInBox; ++iState) {
	    const CFreal currDistance = distances[iState].first;
	    if (MathChecks::isNotZero(currDistance)) {
	      const CFreal weight = 1./currDistance;
	      cf_assert(!MathChecks::isNaN(weight));
	      sumWeights += weight;
	      *states[i] += weight*(*distances[iState].second);

	      counter++;
	    }
	    if (counter == selectionSize) break;
	  }
	  cf_assert(counter > 0);

	  *states[i] /= sumWeights;
	}
	else {
	  // in this case the starting state and the one on which I want to
	  // extrapolate the solution coincide
	  cf_assert(twinState != CFNULL);
	  states[i]->copyData(*twinState);
	}

	countStateIDs++;
	foundBox = true;
	distances.clear();
	cf_assert(distances.capacity() == maxNbStatesInBox);
	break; // get out of the loop over boxes
      }
    }

    cf_assert(foundBox);
  }

  cf_always_assert(countStateIDs == nbStates);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
ShepardMeshInterpolator::needsSockets()
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
