// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PhysicalModel_hh
#define COOLFluiD_Framework_PhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SelfRegistPtr.hh"
#include "PhysicalModelImpl.hh"
#include "Common/CFLog.hh"
#include "Framework/NamespaceStack.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a PhysicalModel.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API PhysicalModel : public Common::NonCopyable<PhysicalModel>,
                      public Config::ConfigObject {

friend class PhysicalModelStack;

public: // methods
 
  /// Set the factory registry
  void setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr);
 
  /// Get the factory registry
  Common::SafePtr<Common::FactoryRegistry> getFactoryRegistry();
  
  /// Set the PhysicalModelImpl once created and configured
  void setPhysicalModelImpl(
  Common::SelfRegistPtr<PhysicalModelImpl> physicalModelImpl);

  /// Set the PhysicalModelImpl
  void unsetPhysicalModelImpl();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args );

  /// Set the total number of equations subsystems
  void setTotalNbEqSS(CFuint totalNbEqSS)
  {
    getEquationSubSysDescriptor().setTotalNbEqSS(totalNbEqSS);
  }

  /// Set the equation variable pattern for each equation subsystem
  void setEquationVarPatterns(const std::vector<std::vector<CFuint> >& eqVarPatterns)
  {
    getEquationSubSysDescriptor().setEquationVarPatterns(eqVarPatterns);
  }

  /// Tells if the equations are solved adimensionalized
  bool isAdimensional() const
  {
    return getImplementor()->isAdimensional();
  }

  /// Set some data describing the current equation subsystem to solve
  /// This is meant to be particularly useful if weak coupling is used
  /// @param iStartVar    first variable ID for the current equation subsystem
  /// @param nbEqsSubSys  number of equations of the current equation subsystem
  /// @param iEqSubSys    ID of the current equation subsystem
  void setEquationSubSysDescriptor(CFuint iStartVar,
      CFuint nbEqsSubSys,
      CFuint iEqSubSys)
  {
    getEquationSubSysDescriptor().set(iStartVar, nbEqsSubSys, iEqSubSys);
  }

  /// Reset to default values some data describing the current equation
  /// subsystem to solve.
  /// This is meant to be particularly useful if weak coupling is used
  void resetEquationSubSysDescriptor()
  {
    // the number of equations is set here to avoid virtual calls
    // in PhysicalModelImpl and allow full inlining of such frequently called
    // function
    getEquationSubSysDescriptor().reset(_nbEquations);
  }

  /// Get some data describing the current equation
  /// subsystem to solve.
  /// This is meant to be particularly useful if weak coupling is used
  EquationSubSysDescriptor& getEquationSubSysDescriptor()
  {
    return getImplementor()->getEquationSubSysDescriptor();
  }

  /// Set current zone name
  void setCurrentZone(const std::string zoneName)
  {
    getImplementor()->setCurrentZone(zoneName);
  }

  /// @todo add method to return a reference
  ///       and an explicit one to return a pointer
  /// @return the PhysicalModelImpl
  Common::SafePtr<PhysicalModelImpl> getImplementor() const
  {
    cf_assert(_setup);
    cf_assert(_physicalModelImpl.isNotNull());
    return _physicalModelImpl.getPtr();
  }

  /// Check if the State is valid
  /// @return true if the State is valid
  bool validate(const State& state)
  {
    cf_assert(_setup);
    cf_assert(_physicalModelImpl.isNotNull());
    return _physicalModelImpl->validate(state);
  }

  /// @return the space dimension of the SubSystem
  CFuint getDim()
  {
    cf_assert(_setup);
    cf_assert(_physicalModelImpl.isNotNull());
    return _dimension;
  }

  /// @return the number of equations of the SubSystem
  CFuint getNbEq()
  {
    cf_assert(_setup);
    cf_assert(_physicalModelImpl.isNotNull());
    return _nbEquations;
  }

  /// @return the name of the PhysicalModel
  std::string getNameImplementor()
  {
    cf_assert(_setup);
    cf_assert(_physicalModelImpl.isNotNull());
    return _nameImplementor;
  }

  /// @return the name of the convective part of the PhysicalModel
  std::string getConvectiveName()
  {
    cf_assert(_setup);
    cf_assert(_physicalModelImpl.isNotNull());
    return _convectiveName;
  }

  /// @return the name of the diffusive part of the PhysicalModel
  std::string getDiffusiveName()
  {
    cf_assert(_setup);
    cf_assert(_physicalModelImpl.isNotNull());
    return _diffusiveName;
  }

  /// @return the name of the linear source part of the PhysicalModel
  std::string getSourceName()
  {
    cf_assert(_setup);
    cf_assert(_physicalModelImpl.isNotNull());
    return _sourceName;
  }

  /// Default destructor
  ~PhysicalModel();

private: // methods

  /// Default constructor
  PhysicalModel(const std::string& name);

private: // member data
  
  /// factory registry to allow polymorphic creation of objects
  Common::SafePtr<Common::FactoryRegistry> m_fr;
  
  /// Physical model implementor
  Common::SelfRegistPtr<PhysicalModelImpl> _physicalModelImpl;

  /// PhysicalModel is setup
  bool _setup;

  /// space dimension of the SubSystem
  CFuint _dimension;

  /// the number of equations of the SubSystem
  CFuint _nbEquations;

  /// the name of this PhysicalModel
  std::string _nameObject;

  /// the name of this PhysicalModel Implementor
  std::string _nameImplementor;

  /// the name of the convective part of this PhysicalModel
  std::string _convectiveName;

  /// the name of the diffusive part of this PhysicalModel
  std::string _diffusiveName;

  /// the name of the linear source part of this PhysicalModel
  std::string _sourceName;

}; // end of class PhysicalModel

//////////////////////////////////////////////////////////////////////////////

class Framework_API PhysicalModelStack : public NamespaceStack<PhysicalModel> {
public:

  /// Returns the instance of the Active PhysicalModel
  /// which is the one on top of the stack
  /// @return SafePtr to the active PhysicalModel
  static Common::SafePtr<PhysicalModel> getActive()
  {
    cf_assert(getInstance().isEnabled());
    return getInstance().top();
  }

  /// Returns the instance of this meshDataStack
  /// This is the access point to the Singleton
  /// @return the instance of the singleton
  static PhysicalModelStack& getInstance()
  {
    static PhysicalModelStack aPhysicalModelStack;
    return aPhysicalModelStack;
  }

protected: // helper functions from NamespaceStack

  /// Gets the name of the PhysicalModel from the Namespace
  /// @param nsp the Namespace from where to get te object name
  std::string getObjectName(const Common::SafePtr<Namespace>& nsp) const;

  /// Creates a PhysicalModel with the supplied name
  /// @param name of the PhysicalModel
  PhysicalModel * createObject(const std::string& name);

}; // end of class PhysicalModelStack;

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PhysicalModel_hh
