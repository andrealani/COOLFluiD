// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_LinearFunctor_hh
#define COOLFluiD_MathTools_LinearFunctor_hh

//////////////////////////////////////////////////////////////////////////////

#include "RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a functor that returns a linear RealVector
/// @author Thomas Wuilbaut
template <CFuint RESULT_SIZE>
class MathTools_API LinearFunctor {
public:

  /// Constructor
  explicit LinearFunctor() :
    _result(RESULT_SIZE)
  {
  }

  /// Default destructor
  ~LinearFunctor()
  {
  }

  void setLinearCoef(const CFreal& value)
  {
    _linearCoef = value;
  }

  void setIndepCoef(const CFreal& value)
  {
    _indepCoef = value;
  }

  /// @return the size of the result
  COOLFluiD::CFuint size() const
  {
    return RESULT_SIZE;
  }

  /// Overloading of operator()
  COOLFluiD::RealVector& operator()(const RealVector& value)
  {
    // must assign it always because some client can have changed it meanwhile
    cf_assert(value.size() == _result.size());
    for(CFuint i=0; i< _result.size(); ++i)
    {
      _result[i] = _indepCoef + (_linearCoef*value[i]);
    }
    CFout << "Result = " << _result[0] << " using value: " << value[0] << "\n";
    return _result;
  }

  private:

  /// the linear coeficient (b) of the equation result= a+bx
  CFreal _linearCoef;

  /// the independent coeficient (a) of the equation result= a+bx
  CFreal _indepCoef;

  /// array storing the temporary solution
  COOLFluiD::RealVector _result;

}; // end class LinearFunctor

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_LinearFunctor_hh
