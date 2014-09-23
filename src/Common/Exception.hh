// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Exception_hh
#define COOLFluiD_Common_Exception_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CodeLocation.hh"
#include "Common/NonCopyable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Manager of behavior of exceptions
class Common_API ExceptionManager : public Common::NonCopyable <ExceptionManager> {
public:

  /// Constructor
  ExceptionManager();
  
  /// Gets the instance of the manager
  static ExceptionManager& getInstance ();

  /// if exception contructor should output
  bool ExceptionOutputs;
  /// if exception contructor should dump backtrace
  bool ExceptionDumps;
  /// if exception contructor should abort execution immedietly 
  bool ExceptionAborts;

}; // class ExceptionManager

//////////////////////////////////////////////////////////////////////////////

/// @brief Base class for all Exceptions in COOLFluiD
///
/// This type extends std::exception, only by providing an implementation for what().
/// <p>what() should optionally describe the situation of the exception in
///   more detail. Interesting information about the exception occurence
///   should not be passed as a part of this string, but rather in separate
///   data members, defined in an appropriate subtype of this class.</p>
/// <p>This class is abstract. Only subclasses can be instantiated.</p>
///
/// @protected
/// <p>Be very carefull with the data members you add in subclasses. There is a great
///   potential for memory leaks. The execution context where the exception
///   is thrown is lost, so if this context holds the only existing reference
///   to objects on the heap, they can no longer be deleted.</p>
///
/// @private
/// <p>The implementation of what() is based on the implementation of
///   std::logic_error. We will not use std::logic_error itself, because the
///   semantics of that type are different, so we do not want our type to be
///   a subtype of std::logic_error.</p>
///
/// @author Tiago Quintino
/// @author Andrea Lani
/// @author Pieter De Ceuninck
/// @author Jan Dockx
class Common_API Exception : public std::exception {
public: // functions

  /// Default copy constructor
  virtual ~Exception () throw ();

  /// Returns a verbose message with all information about this exception
  std::string full_description () const throw ();

  /// Gets the what description string which does not contain EOL's.
  const std::string& str () const throw ();

  /// @return str().c_str();
  const char* what () const throw ();

  /// Append the message to the what() description
  /// @param msg the std::string to be appended
  void append (const std::string& add) throw ();

  /// @returns the Exception name
  /// Pure virtual so must implement this function in each subclass.
  virtual std::string getClassName () throw () { return m_class_name; }

protected: // functions

  /// The constructor is protected to force the developers to create subclasses.
  /// @param msg  A message describing the circumstances of this exception occurence which might be the empty string.
  /// @pre   msg should not contain EOL's
  /// @post  new.what() == what;
  Exception (CodeLocation where, std::string msg, std::string className) throw ();

protected: // data

  /// Similar to m_msg but stores from where the exception was thrown
  /// indicating the file, line and function if possible.
  /// This is an encapsulated data member, which uses value semantics, and
  /// not a reference or a pointer. A reference or a pointer in an
  /// exception leads immediately to a memory leak, so cannot be used.
  CodeLocation m_where;

  /// Stores the message with explanation of what happened
  std::string m_msg;

  /// The subclass name
  std::string m_class_name;

  /// Stores the full description message with explanation of what happened
  /// @post same as full_description()
  std::string m_what;

}; // end class Exception

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_Exception_hh
