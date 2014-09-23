// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_BDF2_InitCN_hh
#define COOLFluiD_Numerics_NewtonMethod_BDF2_InitCN_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a ConvergenceMethod that implements the BDF2_InitCN
/// method for solving non linear systems, where Cranck-Nicholson is
/// used for the first step.
/// @author Thomas Wuilbaut
/// @author Kris Van den Abeele

class BDF2_InitCN : public NewtonIterator {
public:

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit BDF2_InitCN(const std::string& name);

  /// Default destructor
  ~BDF2_InitCN();

  /// Configures the method, by allocating the it's dynamic members.
  /// @param args arguments from where to read the configuration
  virtual void configure ( Config::ConfigArgs& args );

protected: // abstract interface implementations

  /// Take one timestep
  /// @see ConvergenceMethod::takeStep()
  virtual void takeStepImpl();

protected: // variables

  /// A prepare command to use for the first step
  Common::SelfRegistPtr<NewtonIteratorCom> m_prepare1stStep;

  /// An intermediate command to use for the first step
  Common::SelfRegistPtr<NewtonIteratorCom> m_intermediate1stStep;

  /// The string for configuration of m_prepare1stStep command
  std::string m_prepare1stStepStr;

  /// The string for configuration of m_intermediate1stStep command
  std::string m_intermediate1stStepStr;

}; // class BDF2_InitCN

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_BDF2_InitCN_hh
