// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConvectiveVarSet_hh
#define COOLFluiD_Framework_ConvectiveVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Common/OwnedObject.hh"
#include "Common/NonCopyable.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Common/NullableObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a variable set,
/// which provides physical model dependent data and methods
/// associated to a choice of variables
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API ConvectiveVarSet :
      public Common::OwnedObject,
      public Common::NonCopyable<ConvectiveVarSet>,
      public Common::NullableObject
{
public:

  typedef Environment::ConcreteProvider<ConvectiveVarSet,1> PROVIDER;
  typedef Common::SafePtr<Framework::BaseTerm> ARG1;
  friend class Flux;
  
  /// Constructor
  ConvectiveVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /// Default destructor
  virtual ~ConvectiveVarSet();

  /// Set up the private data and give the maximum size of states physical
  /// data to store
  virtual void setup();

  /// Unset up the private data and give the maximum size of states physical
  /// data to store
  virtual void unsetup();
  
  /// Set the PhysicalData corresponding to the given State
  virtual void computePhysicalData (const State& state, RealVector& pdata) = 0;
  
  /// Set the State correspoding to the PhysicalData
  virtual void computeStateFromPhysicalData(const RealVector& pdata, State& state) 
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::computeStateFromPhysicalData()");
  }
  
  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues)
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::computeEigenValues()");
  }
  
  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::getMaxEigenValue()");
  }
  
  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal) 
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::getMaxAbsEigenValue()");
  }
  
  /// Compute the matrix of the right and left eigenvectors and
  /// the vector of the eigenvalues.
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
					 RealMatrix& leftEv, 
					 RealVector& eValues,
					 const RealVector& normal) 
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::computeEigenValuesVectors()");
  }
  
  /// Set jacobian dissipation coefficient
  void setJacobDissipCoeff(CFreal jacobDissip) {_jacobDissip = jacobDissip;}

  /// Set the entity ID to be used in special cases for MHD
  void setEntityID(CFuint entityID)
  {
    _entityID = entityID;
  }
  
  /// Get the size of the extra physical vars
  virtual CFuint getExtraPhysicalVarsSize() {return 0;}
  
  /// Set extra physical variables that will be used before computing the physical data
  /// The default implementation does nothing
  virtual void setExtraPhysicalVars(RealVector* extraVars) {}
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar) 
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::computePerturbedPhysicalData()");
  }
  
  /// Gets the block separator for this variable set
  virtual CFuint getBlockSeparator() const
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::getBlockSeparator()");
    return 0;
  }
  
  /// Gets the average advection vector the specified variable
  /// @param vec the vector to change
  /// @param iVar the variable corresponding to the specified vector
  virtual void getAverageAdvectionVector(RealVector& vec, const CFuint iVar) const
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::getAverageAdvectionVector()");
  }
  
  /// Computes the jacobian matrices
  virtual void computeJacobians()
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::computeJacobians()");
  }
  
  /// Computes the projection of the jacobian matrix to the normal
  virtual void computeProjectedJacobian (const RealVector& normal, RealMatrix& jacob)
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::computeProjectedJacobian()");
  }
  
  /// Computes the scalar part of the jacobian
  virtual void computeScalarJacobian (const RealVector& normal, RealVector& jacob)
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::computeScalarJacobian()");
  }
  
  /// Compute the splitted the jacobian
  virtual void splitJacobian ( RealMatrix& jacobPlus,
                               RealMatrix& jacobMin,
                               RealVector& eValues,
                               const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"ConvectiveVarSet::splitJacobian()");
  }
  
  /// Give dimensional values to the adimensional state variables
  virtual void setDimensionalValues(const State& state, RealVector& result) {result = state;}

  /// Give adimensional values to the dimensional state variables
  virtual void setAdimensionalValues(const State& state, RealVector& result) {result = state;}
  
  /// Set other dimensional values for useful physical quantities
  virtual void setDimensionalValuesPlusExtraValues
  (const State& state, RealVector& result, RealVector& extra)
  {
    setDimensionalValues(state, result);
  }
  
  /// Gets the Class name
  static std::string getClassName() {return "ConvectiveVarSet";}
  
  /// Get variable names
  virtual const std::vector<std::string>& getVarNames() const
  {
    cf_assert(_varNames.size() > 0);
    return _varNames;
  }
  
  /// Get extra variable names
  virtual std::vector<std::string> getExtraVarNames() const
  {
    return std::vector<std::string>();
  }

  /// mask array for convective variables (if = 0, no convection is needed for that variable) 
  const std::vector<bool>& getMaskVariableArray() const
  {
    return _maskArray;
  }

  /// Set the extra data flag
  void setExtraData(bool extraData)
  {
    _extraData = extraData;
  }
  
  /// Tells if there is a source term
  bool hasSourceTerm() const
  {
    return _hasSourceTerm;
  }
  
  /// Get the ID of the current equation subsystem (in weakly coupled simulations)
  CFuint getEqSS() const {return _iEqSubSys;}
  
  /// Set the ID of the current equation subsystem (in weakly coupled simulations)
  void setEqSS(const CFuint iEqSS) {_iEqSubSys = iEqSS;}
  
  /// Set the delta
  void setDelta(CFreal delta) {_delta = delta;}

  /// Returns true if the state doesn't have unphysical values
  /// Default is to do no check and return true.
  /// Derived classes should implement the checks dependent on the concrete physics
  virtual bool isValid (const RealVector& state)  {return true;}
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs) = 0;
  
  /// This nested class is a composition-based adapter functor which provides
  /// an interface to compute the convective physical flux. The actual job is
  /// delegated to the ConvectiveVarSet
  /// @author Andrea Lani
  class Framework_API Flux {

  public: // methods

    /// Constructor
    Flux(ConvectiveVarSet& vs) :
      _vs(vs)
    {
    }

    /// Default destructor
    virtual ~Flux()
    {
    }

    /// Set up the private data
    virtual void setup()
    {
      _vs._fluxArray.resize(PhysicalModelStack::getActive()->getNbEq());
      _vs._physFlux.resize(PhysicalModelStack::getActive()->getNbEq(),PhysicalModelStack::getActive()->getDim());
    }

    /// Unsetup the private data
    virtual void unsetup()
    {
      _vs._fluxArray.resize(0);
      _vs._physFlux.resize(0,0);
    }

    /// @return the size of the result
    CFuint size() const
    {
      return _vs._fluxArray.size();
    }

    /// Overloading of operator()
    RealVector& operator() (const RealVector& pdata,
			    const RealVector& normals)
    {
      _vs.computeFlux(pdata, normals);
      return _vs._fluxArray;
    }
    
    /// Overloading of operator()
    RealMatrix& operator() (const RealVector& pdata)
    {
      _vs.computeStateFlux(pdata);
      return _vs._physFlux;
    }
    
  protected: // data
    
    /// acquaintance of the variable set
    ConvectiveVarSet& _vs;
    
  }; // end class Flux
  
  /// Get the flux
  Flux& getFlux()
  {
    cf_assert(_flux.get() != CFNULL);
    return *_flux;
  }
  
protected: // methods
 
  /// Set the list of the variable names
  void setVarNames ( const std::vector<std::string>& varNames )
  {
    /// @TODO AL check this
    _varNames = varNames;
  }
  
  /// Add a neme to the list of the variable names
  void addVarName(const std::string& varName)
  {
    _varNames.push_back(varName);
  }
  
  /// Set _hasSourceTerm = true
  /// @pre by default _hasSourceTerm = false
  void turnOnSourceTerm ()
  {
    _hasSourceTerm = true;
  }
  
  /// Compute the convective flux, projected in normal
  virtual void computeFlux (const RealVector& pdata, const RealVector& normals) = 0;
  
  /// Compute the physical convective flux
  virtual void computeStateFlux (const RealVector& pdata) = 0;
  
protected: // data
  
  /// number of equations
  CFuint _nbEqs;

  /// the entity ID to be used in special cases related to MHD
  CFuint _entityID;
    
  /// ID of the current equation subsystem (in weakly coupled simulations)
  CFuint _iEqSubSys;
    
  /// extra data flag
  bool _extraData;
  
  /// list of the variables names
  std::vector<std::string> _varNames;
  
  /// mask array for convective variables
  std::vector<bool> _maskArray;
    
  /// this bool is true if there is a source term
  bool _hasSourceTerm;
  
  /// jacobian dissipation coefficient
  CFreal _jacobDissip;

  /// vector to store the temporary positive eigenvalues
  RealVector _eValuesP;
  
  /// vector to store the temporary negative eigenvalues
  RealVector _eValuesM;

  /// array storing the flux
  RealVector _fluxArray;

  /// array storing the physical flux (not projected onto a normal)
  RealMatrix _physFlux;
  
  /// convective flux functor
  std::auto_ptr<Flux> _flux;
  
  /// use for adding jacobian dissipation @see Euler2DCons
  CFreal _delta;
  
}; // end of class ConvectiveVarSet

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(ConvectiveVarSet) // define the factory instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConvectiveVarSet_hh
