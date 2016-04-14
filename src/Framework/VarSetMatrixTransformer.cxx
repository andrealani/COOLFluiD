// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VarSetMatrixTransformer.hh"
#include "PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

VarSetMatrixTransformer::VarSetMatrixTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _isIdentityTransformation(false),
  _transMatrix(),
  _transVec(),
  _transVecVec(),
  _transVecVecMulti()
{
}

//////////////////////////////////////////////////////////////////////////////

VarSetMatrixTransformer::~VarSetMatrixTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void VarSetMatrixTransformer::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);

  _transMatrix.resize(PhysicalModelStack::getActive()->getNbEq(),
		      PhysicalModelStack::getActive()->getNbEq());
  _transVec.resize(PhysicalModelStack::getActive()->getNbEq());
  
  _transVecVecMulti.resize(maxNbTransStates);
  _transVecVec.resize(maxNbTransStates);
  for (CFuint i = 0; i <  _transVecVec.size(); ++i) {
    _transVecVec[i].resize(PhysicalModelStack::getActive()->getNbEq());
    _transVecVec[i] = 0.0;
  }

  _isIdentityTransformation = getIsIdentityTransformation();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _transMatrix = 0.0;
  if (_isIdentityTransformation) {
    // set the transformation matrix to be the IDENTITY matrix
    for (CFuint i = 0; i < nbEqs; ++i) {
      _transMatrix(i,i) = 1.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<State*>* VarSetMatrixTransformer::transform(vector<State*> *const states)
{
  if (!_isIdentityTransformation) {
    const CFuint nbStates = states->size();
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      setMatrix(*(*states)[iState]);
      static_cast<RealVector&>(*_transStateVec[iState]) =
        _transMatrix * (*(*states)[iState]);
    }
    return &_transStateVec;
  }
  return states;
}

//////////////////////////////////////////////////////////////////////////////

State* VarSetMatrixTransformer::transform(State* const state)
{
  if (!_isIdentityTransformation) {
    setMatrix(*state);
    static_cast<RealVector&>(*_transState) = _transMatrix * (*state);
    return _transState;
  }
  return state;
}

//////////////////////////////////////////////////////////////////////////////

vector<RealVector>* VarSetMatrixTransformer::transform(vector<RealVector> *const  states)
{
  if (!_isIdentityTransformation) {
    const CFuint nbStates = states->size();
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      setMatrix((*states)[iState]);
      _transVecVec[iState] = _transMatrix * (*states)[iState];
    }
    return &_transVecVec;
  }
  return states;
}

//////////////////////////////////////////////////////////////////////////////

RealVector* VarSetMatrixTransformer::transform(RealVector* const state)
{
  if (!_isIdentityTransformation) {
    setMatrix(*state);
    _transVec = _transMatrix * (*state);
    return &_transVec;
  }
  return state;
}

//////////////////////////////////////////////////////////////////////////////

vector<State*>* VarSetMatrixTransformer::transformFromRef(vector<State*> *const states)
{
  if (!_isIdentityTransformation) {
    setMatrixFromRef();

    const CFuint nbStates = states->size();
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      static_cast<RealVector&>(*_transStateVec[iState]) =
        _transMatrix * (*(*states)[iState]);
    }
    return &_transStateVec;
  }
  return states;
}

//////////////////////////////////////////////////////////////////////////////

State* VarSetMatrixTransformer::transformFromRef(State* const state)
{
  if (!_isIdentityTransformation) {
    setMatrixFromRef();
    static_cast<RealVector&>(*_transState) = _transMatrix * (*state);
    return _transState;
  }
  return state;
}

//////////////////////////////////////////////////////////////////////////////

vector<RealVector>* VarSetMatrixTransformer::transformFromRef(vector<RealVector> *const states)
{
  if (!_isIdentityTransformation) {
    setMatrixFromRef();
    const CFuint nbStates = states->size();
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _transVecVec[iState] = _transMatrix * (*states)[iState];
    }
    return &_transVecVec;
  }
  return states;
}

//////////////////////////////////////////////////////////////////////////////

RealVector* VarSetMatrixTransformer::transformFromRef(RealVector* const state)
{
  if (!_isIdentityTransformation) {
    setMatrixFromRef();
    _transVec = _transMatrix * (*state);
    return &_transVec;
  }
  return state;
}

//////////////////////////////////////////////////////////////////////////////

vector<RealVector>* VarSetMatrixTransformer::transformMultiFromRef
(vector<RealVector> *const states, CFuint nb)
{
  cf_assert (nb != 0);

  if (!_isIdentityTransformation) {
    setMatrixFromRef();
    const CFuint nbStates = states->size();
    // Eventually resize the transformed vector of RealVector
    for (CFuint i = 0; i < _transVecVecMulti.size(); ++i) {
      _transVecVecMulti[i].resize
	(PhysicalModelStack::getActive()->getNbEq()*nb);
      _transVecVecMulti[i] = 0.0;
    }

    // Transform for
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i<nb; ++i){
      for (CFuint iState = 0; iState < nbStates; ++iState) {
        _transVecVecMulti[iState].slice(i*nb, nbEqs) =
	  _transMatrix * (*states)[iState].slice(i*nb, nbEqs);
      }
    }
    return &_transVecVecMulti;
  }
  return states;
}

//////////////////////////////////////////////////////////////////////////////

std::string VarSetMatrixTransformer::getProviderName
(const std::string& physicalModelName,
 const std::string& first,
 const std::string& second,
 const std::string& third)
{
  if (first == second) return "Identity";
  return (physicalModelName + first + "To" + second + "In" + third);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
