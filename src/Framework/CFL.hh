// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CFL_hh
#define COOLFluiD_Framework_CFL_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Framework/ComputeCFL.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an object where
/// the current CFL value is stored for the current
/// SubSystem.
/// @author Tiago Quintino
/// @author Andrea Lani
class Framework_API CFL :
    public Common::OwnedObject,
    public Config::ConfigObject
{
public:

  /// Default constructor
  CFL();

  /// Default destructor
  virtual ~CFL();

  /// Configure
  virtual void configure ( Config::ConfigArgs& args );

  /// Get the CFL
  CFreal getCFLValue() const { return m_value; }

  /// Set the CFL
  void setCFLValue(const CFreal cfl) { m_value = cfl; }

  /**
   * Update the value of CFL
   * This function takes a convergencestatus pointer. If it is not allocated, 
   * it allocates with subsystem convergence status
   */
  void update(Framework::ConvergenceStatus * cstatus);

  /**
   * Update the value of CFL using the subsystem convergence status
   */
  void update();

  /// Defines the Config Option's of this class
  /// @param opions a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

private:

  /// the CFL value stored
  CFreal m_value;

  /// string for configuration of the CFL term computer
  std::string m_computeCFLStr;

  // CFL calculator
  Common::SelfRegistPtr<ComputeCFL> m_computeCFL;

}; // end of class CFL

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CFL_hh
