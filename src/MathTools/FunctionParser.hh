// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_FunctionParser_hh
#define COOLFluiD_MathTools_FunctionParser_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a parser for analytical functions.
/// Function parser v2.22 by Warp
/// Parses and evaluates the given function with the given variable values.
/// @author  Warp
class MathTools_API FunctionParser {
public:

    /// Constructor
    FunctionParser();

    /// Destructor
    ~FunctionParser();

    /// missing documentation
    /// @param Function missing documentation
    /// @param Vars missing documentation
    /// @return missing documentation
    int Parse(const std::string& Function, const std::string& Vars);

    /// missing documentation
    /// @return missing documentation
    const char* ErrorMsg(void) const;

    /// missing documentation
    /// @param Vars missing documentation
    /// @return missing documentation
    double Eval(const RealVector& Vars);

    /// missing documentation
    /// @return missing documentation
    inline int EvalError(void) const { return EvalErrorType; }

private:

  /// Missing docs
  int varAmount;

  /// Missing docs
  int ParseErrorType;

  /// Missing docs
  int EvalErrorType;

  /// Missing docs
  typedef std::map<std::string, unsigned> VarMap_t;
  VarMap_t Variables;

  /// missing documentation
  struct MathTools_API CompiledCode
  {
    CompiledCode();
    CompiledCode(const CompiledCode&);
    ~CompiledCode();

        unsigned* ByteCode;
        unsigned ByteCodeSize;
        double* Immed;
        unsigned ImmedSize;
        double* Stack;
        unsigned StackSize, StackPtr;
        bool thisIsACopy;
    } Comp;

  /// missing documentation
  /// @return missing documentation
  VarMap_t::const_iterator FindVariable(const char*);

  /// missing documentation
  /// @return missing documentation
  int CheckSyntax(const char*);

  /// missing documentation
  /// @return missing documentation
  bool Compile(const char*);

  /// missing documentation
  /// @return missing documentation
  bool IsVariable(int);

  /// missing documentation
  void AddCompiledByte(unsigned);

  /// missing documentation
  /// @return missing documentation
  int CompileIf(const char*, int);

  /// missing documentation
  /// @return missing documentation
  int CompileElement(const char*, int);

  /// missing documentation
  /// @return missing documentation
  int CompilePow(const char*, int);

  /// missing documentation
  /// @return missing documentation
  int CompileMult(const char*, int);

  /// missing documentation
  /// @return missing documentation
  int CompileAddition(const char*, int);

  /// missing documentation
  /// @return missing documentation
  int CompileComparison(const char*, int);

  /// missing documentation
  /// @return missing documentation
  int CompileAnd(const char*, int);

  /// missing documentation
  /// @return missing documentation
  int CompileOr(const char*, int);

  /// missing documentation
  /// @return missing documentation
  int CompileExpression(const char*, int, bool=false);

  /// missing documentation
  /// @return missing documentation
  FunctionParser(const FunctionParser&);

}; // end class FunctionParser

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_FunctionParser_hh
