// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_VarSetMatrixTransformer_hh
#define COOLFluiD_Framework_VarSetMatrixTransformer_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a matrix variable
/// transformer
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API VarSetMatrixTransformer : public VarSetTransformer {
public:

  typedef Environment::ConcreteProvider<VarSetMatrixTransformer,1> PROVIDER;
  typedef Common::SafePtr<Framework::PhysicalModelImpl> ARG1;

  /// Default constructor without arguments
  VarSetMatrixTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  virtual ~VarSetMatrixTransformer();

  /// Set the data to prepare the simulation
  virtual void setup(const CFuint maxNbTransStates);

  /// Get the transformation matrix
  /// @post the matrix returned is the IDENTITY matrix if no
  ///       transformation has been applied
  const RealMatrix* getMatrix() const
  {
    return &_transMatrix;
  }

  /// Set the transformation matrix
  virtual void setMatrix(const RealVector& state)
  {
  }

  /// Set the transformation matrix from reference values
  virtual void setMatrixFromRef()
  {
    throw Common::NotImplementedException
      (FromHere(), "VarSetMatrixTransformer::setMatrixFromRef()");
  }

  /// Transform a set of state vectors into another one
  std::vector<State*>*
  transform(std::vector<State*> *const states);

  /// Transform a state into another one
  State* transform(State* const state);

  /// Transform a set of state vectors into another one
  std::vector<RealVector>* transform(std::vector<RealVector> *const  states);

  /// Transform a state into another one
  RealVector* transform(RealVector* const state);

  /// Transform a set of states into another set of states
  std::vector<State*>* transformFromRef(std::vector<State*> *const states);

  /// Transform a state into another one
  State* transformFromRef(State* const state);

  /// Transform a set of states into another set of states
  std::vector<RealVector>* transformFromRef
  (std::vector<RealVector> *const states);

  /// Transform a multiple layer of states into another set of states
  std::vector<RealVector>*
  transformMultiFromRef(std::vector<RealVector> *const states, CFuint nb);

  /// Transform a state into another one
  RealVector* transformFromRef(RealVector* const state);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "VarSetMatrixTransformer";
  }

  /// Gets the provider name
  /// @param modelName name of the physical model (convective part)
  /// @param first  first variable set
  /// @param second second variable set
  /// @param third  third variable set
  static std::string getProviderName(const std::string& modelName,
      const std::string& first,
      const std::string& second,
      const std::string& third);

private:

  /// Set the flag telling if the transformation is an identity one
  /// @pre this method must be called during set up
  virtual bool getIsIdentityTransformation() const = 0;
  
  /// Transform a state into another one
  virtual void transform(const Framework::State& state, Framework::State& result)
  {
  }
  
  /// Transform a state into another one from reference precomputed
  /// values (physical data)associated to the given state
  virtual void transformFromRef(const RealVector& pdata, Framework::State& result)
  {
  }

protected:

  /// tells if the transformation is an identity one or not
  bool _isIdentityTransformation;

  /// transformation matrix
  RealMatrix              _transMatrix;

  /// transformed state
  RealVector              _transVec;

  /// transformed vector of states
  std::vector<RealVector> _transVecVec;

  /// transformed vector of states
  std::vector<RealVector> _transVecVecMulti;

}; // end of class VarSetMatrixTransformer

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(VarSetMatrixTransformer)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_VarSetMatrixTransformer_hh
