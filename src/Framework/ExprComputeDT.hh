// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ExprComputeDT_hh
#define COOLFluiD_Framework_ExprComputeDT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ParserException.hh"
#include "MathTools/FunctionParser.hh"
#include "Framework/ComputeDT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class allows the computation of the time step
/// using a vectorial function
/// @author Thomas Wuilbaut

class Framework_API ExprComputeDT : public ComputeDT {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ExprComputeDT(const std::string& name);

  /// Default destructor
  virtual ~ExprComputeDT();

  /// Compute the time step value
  virtual void operator() ();

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

private: // data

  /// Set the Working dir path
  void setFunction();

private:

  /// DT function of iteration
  MathTools::FunctionParser _functionParser;

  /// storage of the residual of the first iteration
  CFreal _firstResidual;

  /// storage of the previous residual
  CFreal _lastResidual;

  /// storage of the maximum residual reached until now
  CFreal _maxResidual;

  /// a temporary vector to store DT
  /// and pass the function parser
  RealVector _eval;

  /// a vector of string to hold the function variables
  std::string _vars;

  /// a string holding the function definition
  std::string _function;

}; // end of class ExprComputeDT

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ExprComputeDT_hh
