// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_VarRegistry_hh
#define COOLFluiD_Framework_VarRegistry_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Framework.hh"
#include "Common/NonCopyable.hh"
#include "Common/GeneralStorage.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class serves as a registry for variables that get dynamically created
/// and therefore cannot be programmed into the code.
/// They are stored as pointers to void, and properly dynamically cast to the
/// correct type upon access.
/// This should be used to store variables with small memory foot print and
/// that are access seldolmly.
/// Each variable is stored with its name given as a std::string and its type,
/// encoded also into std::string.
/// @author Tiago Quintino
class Framework_API VarRegistry : public Common::NonCopyable<VarRegistry> {

public: // functions

  /// Constructor
  VarRegistry();

  /// Destructor
  virtual ~VarRegistry();

  /// regists a variable with the given name
  /// @param name of the variable
  /// @param var pointer to the variable
  template < typename TYPE >
  void registVar ( const std::string& name, TYPE * var )
  {
    if (var == CFNULL)
      throw Common::BadValueException (FromHere(),"Trying to regist null pointer into VarRegistry");
    m_storage.addEntry(name, var );
    m_typestr.addEntry(name, new std::string (DEMANGLED_TYPEID(TYPE)) );
  }

  /// unregists a variable with the given name
  /// @param name of the variable
  template < typename TYPE >
  TYPE * unregistVar ( const std::string& name )
  {
    TYPE * ptr = getVarPtr<TYPE>(name);
    m_storage.removeEntry(name);
    m_typestr.deleteEntry(name);
    return ptr;
  }

  /// accesses the variable
  /// @param name of the variable
  /// @param var pointer to the variable
  template < typename TYPE >
  TYPE& getVar ( const std::string& name )
  {
    TYPE * ptr = getVarPtr<TYPE>(name);
    return *ptr;
  }

  /// accesses the variable
  /// @param name of the variable
  /// @param var pointer to the variable
  template < typename TYPE >
  void setVar ( const std::string& name, const TYPE& value )
  {
    TYPE& var = this->getVar<TYPE>(name);
    var = value;
  }

private: // helper function

  /// accesses the variable
  /// @param name of the variable
  /// @param var pointer to the variable
  template < typename TYPE >
  TYPE * getVarPtr ( const std::string& name )
  {
    void * vptr = m_storage.getEntry(name);
    std::string* n = m_typestr.getEntry(name);
    if (*n != DEMANGLED_TYPEID(TYPE))
    {
      throw Common::BadValueException (FromHere(),"Trying to access variable ["
                               + name + "] of type [" + *n
                               + "] with wrong type ["
                               + DEMANGLED_TYPEID(TYPE) + "]" );
    }
    TYPE * ptr  = static_cast<TYPE*>(vptr);
    return ptr;
  }

private: // data

  /// storage of the variables
  Common::GeneralStorage<void> m_storage;
  /// storage of the variable types
  Common::GeneralStorage<std::string> m_typestr;

}; // end of class VarRegistry

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_VarRegistry_hh
