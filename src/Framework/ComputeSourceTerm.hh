// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeSourceTerm_hh
#define COOLFluiD_Framework_ComputeSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealVector.hh"
#include "Config/ConfigObject.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/MethodStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class ConvectiveVarSet;
    class DiffusiveVarSet;

//////////////////////////////////////////////////////////////////////////////

  /// This class offers an abstract interface to compute
  /// a source term
  /// @author Andrea Lani
  /// @author Tiago Quintino
template < typename METHODDATA >
class ComputeSourceTerm : public MethodStrategy<METHODDATA> {

public:

  typedef Framework::BaseMethodStrategyProvider
  <METHODDATA,ComputeSourceTerm<METHODDATA> > PROVIDER;

  /// Constructor
  ComputeSourceTerm(const std::string& name) :
    MethodStrategy<METHODDATA>(name),
    _sockets(),
    _globalSockets(),
    _iVar(0),
    _useAnalyticalJacob(false),
    _analyticalJacob(false),
    _isPerturb(false)
  {
    this->addConfigOptionsTo(this);
    this->setParameter("UseAnalyticalJacob",&_useAnalyticalJacob);
  }

  /// Default destructor
  virtual ~ComputeSourceTerm()
  {
  }

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options)
  {
    options.template addConfigOption< bool >
      ("UseAnalyticalJacob",
      "Flag forcing to use the analytical jacobian");
  }

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup() 
  {
    MethodStrategy<METHODDATA>::setup();
  } 
  
  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void unsetup() 
  {
    MethodStrategy<METHODDATA>::unsetup();
  }
  
  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Config::ConfigObject::configure(args);
  }
  
  /// Set the flag to tell if the source term is being perturbed
  void setPerturb(bool isPerturb)
  {
    _isPerturb = isPerturb;
  }
  
  /// Compute the source term
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian) = 0;
  
  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<BaseDataSocketSink> > sockets = _sockets.getAllSinkSockets();
    std::vector<Common::SafePtr<BaseDataSocketSink> > globalSockets = _globalSockets.getAllSinkSockets();

    for (CFuint i = 0; i < globalSockets.size(); ++i) {
      sockets.push_back(globalSockets[i]);
    }

    return sockets;
  }
  
  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > sockets = _sockets.getAllSourceSockets();
    return sockets;
  }
  
  /// Tells if the source term has an analytical jacobian
  bool useAnalyticalJacob() const
  {
    return _useAnalyticalJacob;
  }

  /// Set the analtyical jacobian matrix
  void setAnalyticalJacob(bool flag)
  {
    _analyticalJacob = flag;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeSourceTerm";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
protected:

  /// Tells if the source term is being perturbed
  bool isPerturb() const
  {
    return _isPerturb;
  }

protected: //data

  /// The set of data sockets to be used by the strategy
  Framework::DynamicDataSocketSet<> _sockets;
  /// The set of data sockets to be used by the strategy
  Framework::DynamicDataSocketSet<Framework::GLOBAL> _globalSockets;
  
  /// which variable is currently being perturbed
  CFuint _iVar;
  
  /// flag telling if the analytical jacobian has to be used
  bool _useAnalyticalJacob;
  
  /// flag telling if the analytical jacobian has to be computed
  bool _analyticalJacob;
  
  /// flag telling if the source term is being perturbed
  bool _isPerturb;
  
}; // end of class ComputeSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeSourceTerm_hh
