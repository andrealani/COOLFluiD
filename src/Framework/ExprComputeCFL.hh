// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ExprComputeCFL_hh
#define COOLFluiD_Framework_ExprComputeCFL_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeCFL.hh"
#include "MathTools/FunctionParser.hh"
#include "Common/ParserException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class allows the computation of the CFL using a vectorial function
/// @author Tiago Quintino
/// @author Andrea Lani
class Framework_API ExprComputeCFL : public ComputeCFL {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ExprComputeCFL(const std::string& name);

  /// Default destructor
  ~ExprComputeCFL();

  /// Check if the stop condition has been achieved
  void operator() (const ConvergenceStatus& cstatus);

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

private: // methods

  /// Setup the FunctionParser that will parse the expression for CFL.
  void setFunction();

private: // data

  /// CFL function of iteration
  MathTools::FunctionParser _functionParser;

  /// storage of the residual of the first iteration
  CFreal _firstResidual;

  /// storage of the previous residual
  CFreal _lastResidual;

  /// storage of the maximum residual reached until now
  CFreal _maxResidual;

  /// a temporary vector to store CFL
  /// and pass the function parser
  RealVector _eval;

  /// a vector of string to hold the functions
  std::string _vars;

  /// a string holding the function definition
  std::string _function;

}; // end of class ExprComputeCFL

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ExprComputeCFL_hh
