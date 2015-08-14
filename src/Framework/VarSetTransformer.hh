// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_VarSetTransformer_hh
#define COOLFluiD_Framework_VarSetTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"
#include "Common/NullableObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Common/NotImplementedException.hh"
#include "Common/OwnedObject.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a vector variable
/// transformer
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API VarSetTransformer : public Common::OwnedObject,
					public Common::NullableObject {
 public:
  
  typedef Environment::ConcreteProvider<VarSetTransformer, 1> PROVIDER;
  typedef Common::SafePtr<Framework::PhysicalModelImpl> ARG1;
  
  /// Default constructor without arguments
  VarSetTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /// Default destructor
  virtual ~VarSetTransformer();
  
  /// Set the data to prepare the simulation
  virtual void setup(const CFuint maxNbTransStates);
  
  /// Transform an array of variables into another one of potentially different length
  virtual void transform(const RealVector& state, RealVector& result)
  {
    throw Common::NotImplementedException(FromHere(), "VarSetTransformer::transform()");
  }
  
  /// Transform a set of state vectors into another one
  virtual std::vector<State*>* transform(std::vector<State*> *const  states)
  {
    const CFuint nbStates = states->size();
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _iState = iState;
      transform((const Framework::State&)*(*states)[iState], 
		(Framework::State&)*_transStateVec[iState]);
    }
    return &_transStateVec;
  }
  
  /// Transform a state into another one
  virtual State* transform(State* const state)
  {
    transform((const Framework::State&)*state, (Framework::State&)*_transState);
    return _transState;
  }
  
  /// Transform a state into another one from reference precomputed
  /// values (physical data)associated to the given state
  virtual State* transformFromRefData(const RealVector& pdata)
  {
    transformFromRef(pdata, *_transState);
    return _transState;
  }
  
  /// Transform a set of state vectors into another one from reference precomputed
  /// values (physical data)associated to the given state
  virtual std::vector<State*>* transformFromRefData(std::vector<RealVector> *const pdata)
  {
    const CFuint nbStates = pdata->size();
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _iState = iState;
      transformFromRef((*pdata)[iState], *_transStateVec[iState]);
    }
    return &_transStateVec;
  }
  
  /// Set a pointer to an array of state vectors with extra values
  /// @pre if != CFNULL extra values will be computed per degree of freedom
  virtual void setExtraValues(std::vector<RealVector> *const extraValues)
  {
    _extraValues = extraValues;
  }
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "VarSetTransformer ";
  }

  /// Gets the provider name
  /// @param modelName name of the physical model (convective part)
  /// @param first     first variable set
  /// @param second    second variable set
  static std::string getProviderName(const std::string& modelName,
				     const std::string& first,
				     const std::string& second);
  
  /// Set the local state ID
  void setLocalID(CFuint localStateID) {_localStateID = localStateID;}
   
protected:
  
  /// Transform a state into another one
  virtual void transform(const Framework::State& state, Framework::State& result) = 0;
  
  /// Transform a state into another one from reference precomputed
  /// values (physical data)associated to the given state
  virtual void transformFromRef(const RealVector& pdata, Framework::State& result) = 0;
  
  /// Get the local state ID
  CFuint getLocalID() const {return _localStateID;}
  
protected: // data

  /// local ID of the state
  CFuint _iState;

  /// local ID in the whole mesh
  CFuint _localStateID;

  /// arrays of extra values
  std::vector<RealVector>* _extraValues;

  /// transformed state
  State*               _transState;

  /// transformed vector of states
  std::vector<State*> _transStateVec;

}; // end of class VarSetTransformer

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(VarSetTransformer)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_VarSetTransformer_hh
