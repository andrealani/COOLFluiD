// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DiffusiveVarSet_hh
#define COOLFluiD_Framework_DiffusiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Common/NullableObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Common/OwnedObject.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a diffusive
/// variable set, which provides physical model dependent data
/// and methods associated to a choice of variables
/// @author Andrea Lani
class Framework_API DiffusiveVarSet : public Common::OwnedObject,
                        public Config::ConfigObject,
                        public Common::NullableObject {
public:

  typedef Environment::ConcreteProvider<DiffusiveVarSet,2> PROVIDER;
  typedef const std::string& ARG1;
  typedef Common::SafePtr<Framework::PhysicalModelImpl> ARG2;

  /// Constructor
  DiffusiveVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  virtual ~DiffusiveVarSet();

  /// Get the diffusive flux
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius)
  {
    throw Common::NotImplementedException (FromHere(),"DiffusiveVarSet::setFlux()");
    return _flux;
  }

  /// Get the diffusive flux
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius)
  {
    throw Common::NotImplementedException (FromHere(),"DiffusiveVarSet::setFlux()");
    return _fluxVec;
  }
  
  /// Get the diffusive flux
  virtual void computeFluxJacobian(const RealVector& state,
				   const RealVector& gradientJacob,
				   const RealVector& normal,
				   const CFreal& radius,
				   RealMatrix& fluxJacob)
  {
    throw Common::NotImplementedException (FromHere(),"DiffusiveVarSet::computeFluxJacobian()");
  }
  
  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    CFAUTOTRACE;

    ConfigObject::configure(args);

    // add here configuration, specific of this class
  }

  /// set up
  virtual void setup()
  {
    CFAUTOTRACE;
    _flux.resize(PhysicalModelStack::getActive()->getNbEq());
    _flux = 0.;
    _fluxVec.resize(PhysicalModelStack::getActive()->getNbEq(),
		    PhysicalModelStack::getActive()->getDim());
    _fluxVec = 0.;
  }

  /// Compute physical data associated to the given states
  virtual void computeStatesData(const std::vector<State*>& states,
				 const CFuint nbStatesInVec)
  {
    throw Common::NotImplementedException (FromHere(),"VarSet::computeStatesData()");
  }
  
  /// Compute the left perturbed states data
  virtual void computePerturbedStatesData
  (const std::vector<Framework::State*>& states,
  const CFuint nbStatesInVec,
  const CFuint iVar)
  {
    throw Common::NotImplementedException (FromHere(),"VarSet::computePerturbedLeftStatesData()");
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "DiffusiveVarSet";
  }

  /// Get variable names
  const std::vector<std::string>& getVarNames() const
  {
    cf_assert(_varNames.size() > 0);
    return _varNames;
  }

  /// Set the flag forcing to freeze the diffusive coefficients
  void setFreezeCoeff(bool freezeCoeff)
  {
    _freezeDiffCoeff = freezeCoeff;
  }

  /// Get the flag forcing to freeze the diffusive coefficients
  bool getFreezeCoeff() const
  {
    return _freezeDiffCoeff;
  }
  
  /// set the coordinate of the current face to be processed
  /// @pre faceCoord must be used as an array with size equal to the space dimension
  void setFluxCoord(CFreal *const faceCoord)
  {
    _faceCoord = faceCoord; 
  }
  
  /// Set the entity ID to be used in special cases for MHD
  void setEntityID(CFuint entityID)
  {
    _entityID = entityID;
  }
  
protected: // methods

  /// Set the list of the variable names
  void setVarNames(const std::vector<std::string>& varNames)
  {
    if (_varNames.size() == 0) {
      _varNames = varNames;
    }
  }

  /// Add a name to the list of the variable names
  void addVarName(const std::string& varName)
  {
    _varNames.push_back(varName);
  }

private: // methods

  /// Private Copy Constructor
  DiffusiveVarSet(const DiffusiveVarSet& v);

  /// Private Assignement operator
  const DiffusiveVarSet& operator=(const DiffusiveVarSet& v);

protected:

  /// the entity ID to be used in special cases related to MHD
  CFuint _entityID;
  
  /// list of the variables names
  std::vector<std::string> _varNames;

  /// space coordinates of the face center
  CFreal* _faceCoord;
  
  /// flux
  RealVector _flux;

  /// flux vector
  RealMatrix _fluxVec;

  /// Various Indexes
  CFuint _iState;
  CFuint _jState;

  /// flag telling to freeze the diffusive coefficients
  bool _freezeDiffCoeff;

}; // end of class DiffusiveVarSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(DiffusiveVarSet) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DiffusiveVarSet_hh
