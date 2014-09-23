// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SERComputeCFL_hh
#define COOLFluiD_Framework_SERComputeCFL_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeCFL.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a switched evolution relaxation algorithm   
/// @author Khalil Bensassi
class Framework_API SERComputeCFL : public ComputeCFL {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  SERComputeCFL(const std::string& name);

  /// Default destructor
  ~SERComputeCFL();

  /// Check if the stop condition has been achieved
  void operator() (const ConvergenceStatus& m_cstatus);

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

private:
  /// minimum Value of the CFL
  CFreal  m_Min;
  /// Initial Value of the CFL
  CFreal m_initialCFL;
  /// storage of the residual on the second iteration
  CFreal  _Res0;
  /// coefficient for the dynamic CFL calculation
  CFreal m_coeffCFL;
 /// Absolute Minimum Values 
  CFreal m_maxCFL;
  /// power value
  CFreal m_powerExp;
  /// maximum allowable CFL number
  CFreal m_LimitCFL;
  /// bool variable that  allow to switch between local 
  /// and global minimum creteria 
  bool m_Tol;



}; // end of class SERComputeCFL

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SERComputeCFL_hh
