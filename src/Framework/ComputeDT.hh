// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeDT_hh
#define COOLFluiD_Framework_ComputeDT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class defines an object to compute the time step DT
/// @author Thomas Wuilbaut
class Framework_API ComputeDT : public Common::OwnedObject,
                                public Config::ConfigObject {
public:

  typedef Environment::ConcreteProvider<ComputeDT,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Constructor
  ComputeDT(const std::string& name);

  /// Default destructor
  virtual ~ComputeDT();

  /// Compute the time step value
  virtual void operator() () = 0;

  /// Configure
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName() { return "ComputeDT"; }

}; // end of class ComputeDT

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(ComputeDT)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeDT_hh
