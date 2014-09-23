// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_VectorialFunction_hh
#define COOLFluiD_Framework_VectorialFunction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ParserException.hh"
#include "Common/StringOps.hh"
#include "MathTools/RealVector.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools { class FunctionParser; }

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Function that defines the values
/// for a vector field.
/// @author Tiago Quintino
class Framework_API VectorialFunction {

public: // functions

  /// Default constructor without arguments
  VectorialFunction();

  /// Default destructor
  ~VectorialFunction();

  /// Evaluate the Vectorial Function given the values of the variables.
  /// @param vars values of the variables to substitute in the function.
  /// @param value the placeholder vector for the result
  void evaluate(const RealVector& varValue, RealVector& value) const;

  /// Evaluate the Vectorial Function given the values of the variables.
  /// @param iVar index of the variable to compute
  /// @param vars values of the variables to substitute in the function.
  /// @param value the result
  void evaluate(CFuint iVar, const RealVector& varValue, CFreal& value);

  /// Evaluate the Vectorial Function given the values of the variables
  /// and return it in the stored result. This function allows this class to work
  /// as a functor.
  /// @param vars values of the variables to substitute in the function.
  RealVector& operator()(const RealVector& varValue);

  /// @return if the VectorialFunctionParser has been parsed yet.
  bool isParsed() const
  {
    return m_isParsed;
  }

  /// Set Function
  void setFunctions(const std::vector<std::string>& functions);

  /// Set Vars
  void setVariables(const std::vector<std::string>& vars);

  /// Parse the strings to extract the functions for each line of the vector.
  /// @param functions vector of string describing the functions for each line.
  /// @param vars the variables to be taken into account.
  /// @throw ParserException if there is an error while parsing
  void parse();

  /// Gets the number of variables
  /// @pre only call this function if already parsed
  /// @returns the number of varibles
  CFuint getNbVars() const
  {
    cf_assert ( isParsed() );
    return m_nbVars;
  }

  /// Gets the number of variables
  /// @pre only call this function if already parsed
  /// @returns the number of varibles
  CFuint getNbFuncs() const
  {
    cf_assert ( isParsed() );
    return m_functions.size();
  }

protected: // helper functions

  /// Clears the m_parsers deallocating the memory.
  void clear();

private: // data

  /// flag to indicate if the functions have been parsed
  bool m_isParsed;

  /// vector holding the names of the variables
  std::string m_vars;

  /// number of variables
  CFuint m_nbVars;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// vector holding the parsers, one for each entry in the vector
  std::vector<MathTools::FunctionParser*> m_parsers;

  /// storage of the result for using the class as functor
  RealVector m_result;

}; // end of class VectorialFunction

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_VectorialFunction_hh
