// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeCFL_hh
#define COOLFluiD_Framework_ComputeCFL_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"
#include "Framework/ConvergenceStatus.hh"
#include "Framework/SubSystemStatus.hh"
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class CFL;

//////////////////////////////////////////////////////////////////////////////

/// This class defines the condition to be satisfied to stop
/// the computation
/// @author Andrea Lani
class Framework_API ComputeCFL : public Common::OwnedObject,
                   public Config::ConfigObject {
public:

  typedef Environment::ConcreteProvider<ComputeCFL,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Constructor
  ComputeCFL(const std::string& name) :
    Common::OwnedObject(),
    ConfigObject(name)
  {
  }

  /// Default destructor
  virtual ~ComputeCFL()
  {
  }

  /// Check if the stop condition has been achieved
  virtual void operator() (const Framework::ConvergenceStatus& m_cstatus) = 0;

  /// Configure
  virtual void configure ( Config::ConfigArgs& args )
  {
    ConfigObject::configure(args);
  }

  /// Set the CFL
  void setCFL(Common::SafePtr<Framework::CFL> cfl)
  {
    _cfl = cfl;
  }

  /// Gets the Class name
  static std::string getClassName() {   return "ComputeCFL";  }

protected:

  Common::SafePtr<Framework::CFL> _cfl;

}; // end of class ComputeCFL

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(ComputeCFL) // declare the factory of this type

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeCFL_hh
