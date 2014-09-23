// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_InteractiveComputeDT_hh
#define COOLFluiD_Framework_InteractiveComputeDT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"
#include "ComputeDT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class allows to change the DT interactively
/// @author Thomas Wuilbaut
/// @author Tiago Quintino
class Framework_API InteractiveComputeDT : public ComputeDT {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  InteractiveComputeDT(const std::string& name);

  /// Default destructor
  virtual ~InteractiveComputeDT();

  /// Computes the time step value
  void operator() ();

  /// Configures the options for this object.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure ( Config::ConfigArgs& args );

private:

  /// current interactive value of the DT
  CFreal m_dt_value;

}; // end of class InteractiveComputeDT

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_InteractiveComputeDT_hh
