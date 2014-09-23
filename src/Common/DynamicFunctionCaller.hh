// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_DynamicFunctionCaller_hh
#define COOLFluiD_Common_DynamicFunctionCaller_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown
/// when a dynamic function is called and fails
/// @author Tiago Quintino
class Common_API DynamicFunctionException : public Common::Exception {
public:

  /// Constructor
  DynamicFunctionException (const Common::CodeLocation& where, const std::string& what);

}; // end of class DynamicFunctionException

//////////////////////////////////////////////////////////////////////////////

/// This class allows to call void functions on classes that derive from it,
/// dynamically by simply passing the name of the function to execute.
/// @author Tiago Quintino
template < typename TYPE >
class DynamicFunctionCaller
{
  public:

  virtual ~DynamicFunctionCaller(){}

  protected:

  /// pointer to dynamic function
  typedef void (TYPE::*dynamic_function)();

  /// type that maps a std::string with a function
  typedef std::map<std::string,dynamic_function,std::less<std::string> > MapString2Func;

  protected:

    /// Gets the pointer to the function
    /// @post returns CFNULL is not found
    dynamic_function get_dynamic_function(const std::string & func)
    {
      typename MapString2Func::iterator key_pair = db.find(func);
      dynamic_function dfunc = CFNULL;
      if (key_pair != db.end()) { dfunc = key_pair->second; }
      return dfunc;
    }

    /// Runs the function
    /// @post returns CFNULL if no function was found
    dynamic_function run_dynamic_function(const std::string & func)
    {
      dynamic_function dfunc = get_dynamic_function(func);
      if ( dfunc != CFNULL )
      {
         TYPE * ptr = dynamic_cast<TYPE*>(this);
         (ptr->*dfunc)();
      }
      return dfunc;
    }

    /// Adds the function to the database
    void add_dynamic_function(const std::string & func, dynamic_function fp)
    {
      typename MapString2Func::iterator key_pair = db.find(func);

      // check if function doesnt already exist if yes add it
      if (key_pair == db.end())
      {
        db[func] = fp;
      }
      else // throw an exception
      {
        std::string msg("function already added: ");
        msg += func;
        throw DynamicFunctionException (FromHere(),msg);
      }
    }

  /// Executes function on class or on parent if not found
  template <typename PARENT>
  void hierarchical_run_function(const std::string & func)
  {
    dynamic_function dfunc = DynamicFunctionCaller<TYPE>::run_dynamic_function(func);
    if ( dfunc == CFNULL )
    {
        PARENT* obj = dynamic_cast<PARENT*>(this);
        obj->PARENT::run_function(func);
    }
  }

  /// Function that declares which functions can be called
  virtual void build_dynamic_functions() = 0;

  private:
    /// database of dynamic functions to be run
    MapString2Func db;
};

//////////////////////////////////////////////////////////////////////////////

} // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Common_DynamicFunctionCaller_hh
