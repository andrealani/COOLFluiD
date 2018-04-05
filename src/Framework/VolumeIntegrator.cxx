// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeIntegrator.hh"
#include "PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

VolumeIntegrator::VolumeIntegrator() :
  Integrator<VolumeIntegratorImpl>(),
  _intgGeo(CFNULL),
  _intgSol(CFNULL),
  _solValues(),
  _coord(),
  _gradValues(),
  _jacob(),
  _detJacobian()
{
}

//////////////////////////////////////////////////////////////////////////////

VolumeIntegrator::~VolumeIntegrator()
{
  for(CFuint i = 0; i < _solValues.size(); ++i) {
    deletePtr(_solValues[i]);
  }

  for(CFuint i = 0; i < _coord.size(); ++i) {
    deletePtr(_coord[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolumeIntegrator::buildCoordinatesAndStates(const IntegratorPattern& pattern)
{
  const CFuint maxSize = pattern.totalNbPts();

  // clear and resize the storage for temporary interpolated solution values
  for(CFuint i = 0; i < _solValues.size(); ++i) {
    deletePtr(_solValues[i]);
  }
  vector<State*>().swap(_solValues);

  _solValues.resize(maxSize);

  for(CFuint i = 0; i < _solValues.size(); ++i) {
    _solValues[i] = new State();
  }

  // clear and resize the storage for temporary interpolated coordinates
  const bool notOwnedByState = true;
  for(CFuint i = 0; i < _coord.size(); ++i) {
    deletePtr(_coord[i]);
  }
  vector<Node*>().swap(_coord);

  _coord.resize(maxSize);
  for(CFuint i = 0; i < _coord.size(); ++i) {
    _coord[i] = new Node(notOwnedByState);
  }

  for(CFuint i = 0; i < maxSize; ++i){
    _solValues[i]->setSpaceCoordinates(_coord[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void VolumeIntegrator::setup()
{
  CFAUTOTRACE;
  
  Integrator<VolumeIntegratorImpl>::setup();
  
  const IntegratorPattern pattern = getMaxIntegratorPattern();

  CFLogDebugMin("Maximum integrator pattern: " << pattern << "\n");

  buildCoordinatesAndStates(pattern);

  const CFuint maxSize = pattern.totalNbPts();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  // resize the storage for temporary interpolated solution gradient values
  _gradValues.resize(maxSize);
  for(CFuint i = 0; i < maxSize; ++i){
    _gradValues[i].resize(pattern.getNbShapeFunctions(),nbDim);
  }

  // resize the storage for temporary jacobian
  _jacob.resize(maxSize);
  for(CFuint i = 0; i < maxSize; ++i){
    _jacob[i].resize(nbDim,nbDim);
  }

  // resize the storage for temporary determinants of jacobians
  _detJacobian.resize(maxSize);
  _detJacobian = 0.;
}

//////////////////////////////////////////////////////////////////////////////

void VolumeIntegrator::setNbSolQuadraturePoints(GeometricEntity *const geo,
						CFuint& nbPoints)
{
  nbPoints = getSolutionIntegrator(geo)->getCoeff().size();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
