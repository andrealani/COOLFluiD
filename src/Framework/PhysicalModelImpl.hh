// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PhysicalModelImpl_hh
#define COOLFluiD_Framework_PhysicalModelImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Common/NonCopyable.hh"

#include "Common/OwnedObject.hh"

#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"

#include "Config/ConfigObject.hh"

#include "Environment/ConcreteProvider.hh"

#include "Framework/EquationSubSysDescriptor.hh"

/// @todo broken after release 2009.3 (added to fix State class missing error)
#include "Framework/State.hh" 

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
/// @todo broken after release 2009.3
    class State;
    class BaseTerm;
    class PhysicalPropertyLibrary;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a PhysicalModelImpl.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API PhysicalModelImpl : public Common::OwnedObject,
                          public Config::ConfigObject,
                          public Common::NonCopyable<PhysicalModelImpl> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  typedef Environment::ConcreteProvider<PhysicalModelImpl,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Constructor without arguments
  PhysicalModelImpl(const std::string& name);

  /// Default destructor
  virtual ~PhysicalModelImpl();

  /// Set up
  virtual void setup();

  /// Get the convective name
  virtual std::string getConvectiveName() const = 0;

  /// Get the diffusive name
  virtual std::string getDiffusiveName() const = 0;

  /// Get the Source name
  virtual std::string getSourceName() const { return "Null"; }

  /// @return the space dimension of the SubSystem
  virtual CFuint getDimension() const = 0;

  /// @return the number of equations of the SubSystem
  virtual CFuint getNbEquations() const = 0;
  
  /// @return the number of equations of the SubSystem
  bool is2DHalf() const {return m_is2DHalf;}
  
  /// Check if this state is in a valid state
  virtual bool validate(const State& state) const  {  return true; }

  /// Get the convective term
  virtual Common::SafePtr<BaseTerm> getConvectiveTerm() const = 0;

  /// Get the diffusive term
  virtual Common::SafePtr<BaseTerm> getDiffusiveTerm() const = 0;

  /// Get the source term
  virtual Common::SafePtr<BaseTerm> getSourceTerm() const = 0;

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()  {  return "PhysicalModelImpl"; }

  /// Get the jacobians (accessor)
  const std::vector<RealMatrix>& getJacobians() const  {  return _jacobians; }

  /// Get the jacobians (mutator)
  std::vector<RealMatrix>* getJacobians()  {   return &_jacobians;  }

  /// Get the reference values to adimensionalize the state vector
  const RealVector& getRefStateValues() const {   return _refStateValues;  }

  /// Get the reference length for geometric scaling
  CFreal getRefLength() const  {    return _refLength;  }

  /// Get the reference time for time scaling
  CFreal getRefTime() const  {    return _refTime;   }

  /// Tells if the equations are solved adimensionalized
  bool isAdimensional() const  {    return _isAdimensional;   }

  /// Get some data describing the current equation
  /// subsystem to solve.
  /// This is meant to be particularly useful if weak coupling is used
  EquationSubSysDescriptor& getEquationSubSysDescriptor() {  return _eqSubSysDescriptor; }

  /// @return the library computing the physical properties
  /// @post the pointer will be casted to the given library
  ///       type pointer
  template <class LIBRARY>
  Common::SafePtr<LIBRARY> getPhysicalPropertyLibrary() const
  {
    cf_assert(_physicalPropLib.isNotNull());
    return dynamic_cast<LIBRARY*>(_physicalPropLib.getPtr());
  }

  /// Set current zone name
  virtual void setCurrentZone(const std::string zoneName)  {  _currentZoneName = zoneName;  }

private: // methods

  /// Private Copy Constructor
  PhysicalModelImpl(const PhysicalModelImpl& v);

  /// Private Assignement operator
  const PhysicalModelImpl& operator=(const PhysicalModelImpl& v);

  /// Set the reference values
  /// Default empty implementation
  virtual void setReferenceValues()
  {
  }

  /// Set the reference value for time
  /// By default set to the same value as the reference length
  /// For convection equations --> should still be okay this way if only reference length is set.
  virtual void setReferenceTime()  {  _refTime = getRefLength(); }

  /// Set the physical data,
  /// Default empty implementation
  virtual void computePhysicalData() {}

protected: // data

  /// jacobians of this physical model (one for each space dimension)
  std::vector<RealMatrix> _jacobians;

  /// flag telling if the equations are solved adimensional
  bool _isAdimensional;

  /// for simulations for which the properties change from one zone to the other:
  /// string containing the name of the current zone
  std::string _currentZoneName;

  /// reference time for time scaling
  CFreal _refTime;

private: // data

  /// vector of the reference values (to scale the state vector)
  /// used at configuration time
  std::vector<CFreal> _refStateValuesConf;

  /// set of the reference values (to scale the state
  /// vector when needed) used for computation
  RealVector _refStateValues;

  /// reference length for geometric scaling
  CFreal _refLength;

  /// descriptor of the current subsystem of equations
  EquationSubSysDescriptor _eqSubSysDescriptor;

  /// physical properties library
  Common::SelfRegistPtr<Framework::PhysicalPropertyLibrary> _physicalPropLib;

  /// name of the physical properties library
  std::string _physicalPropLibStr;
  
  /// flag for telling that the model is 2D and 1/2
  bool m_is2DHalf;
  
}; // end of class PhysicalModelImpl

//////////////////////////////////////////////////////////////////////////////

  } // namespace PhysicalModelImpl
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(PhysicalModelImpl)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PhysicalModelImpl_hh
