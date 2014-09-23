// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ContourIntegrator.hh"
#include "PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

ContourIntegrator::ContourIntegrator() :
  Integrator<ContourIntegratorImpl>(),
  _unitFaceNormals(CFNULL),
  _faceJacobian(),
  _values(),
  _coord()
{
}

//////////////////////////////////////////////////////////////////////////////

ContourIntegrator::~ContourIntegrator()
{
  for(CFuint i = 0; i < _values.size(); ++i) {
    deletePtr(_values[i]);
  }

  for(CFuint i = 0; i < _coord.size(); ++i) {
    deletePtr(_coord[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ContourIntegrator::buildCoordinatesAndStates
(const IntegratorPattern& pattern)
{
  const CFuint maxSize = pattern.totalNbPts();

  // clear and resize the storage for temporary interpolated solution values
  for(CFuint i = 0; i < _values.size(); ++i) {
    deletePtr(_values[i]);
  }
  vector<State*>().swap(_values);

  _values.resize(maxSize);
  for(CFuint i = 0; i < _values.size(); ++i) {
    _values[i] = new State();
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
    _values[i]->setSpaceCoordinates(_coord[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ContourIntegrator::setup()
{
  CFAUTOTRACE;

  Integrator<ContourIntegratorImpl>::setup();

  const IntegratorPattern pattern = getMaxIntegratorPattern();

  CFLogDebugMin("Maximum integrator pattern: " << pattern << "\n");

  buildCoordinatesAndStates(pattern);

  // resize the storage for temporary face determinants
  _faceJacobian.resize(pattern.nbSteps());
  for (CFuint i = 0; i < pattern.nbSteps(); ++i) {
    _faceJacobian[i].resize(pattern.nbPts(i));
    _faceJacobian[i] = 0.;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ContourIntegrator::setNbSolQuadraturePoints
(const GeometricEntity* geo, CFuint& nbPoints)
{
  const vector<RealVector>& coeff = getSolutionIntegrator(geo)->getCoeff();

  CFuint counter = 0;
  for (CFuint ic = 0; ic < coeff.size(); ++ic)
  {
    counter += coeff[ic].size();
  }

  nbPoints = counter;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
