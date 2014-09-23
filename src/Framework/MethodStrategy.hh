// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MethodStrategy_hh
#define COOLFluiD_Framework_MethodStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SharedPtr.hh"
#include "Framework/NumericalStrategy.hh"
#include "Common/NullableObject.hh"
#include "Framework/BaseMethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

  /// This class represent a Strategy to accomplish some action in a Method.
  /// It is a template intended to function as the command that it is
  /// parametrized with the DATA object shared between all MethodStrategy's
  /// belonging to the same Method.
  /// @author Tiago Quintino
template <class DATA>
class MethodStrategy : public NumericalStrategy,
                       public Common::NullableObject {
public:

  /// Constructor.
  explicit MethodStrategy(const std::string& name) :
    NumericalStrategy(name),
    Common::NullableObject()
  {
  }

  /// Default destructor
  virtual ~MethodStrategy() {}

  /// Gets the shared DATA object for the Method.
  /// @return reference to the Method's Data
  DATA& getMethodData()
  {
    cf_assert(m_dataPtr.isNotNull());
    return *m_dataPtr;
  }

  /// Sets the pointer to the shared DATA object
  /// @param dataPtr is the pointer to the Data to be set.
  void setMethodData(const Common::SharedPtr<DATA>& dataPtr)  
  {  m_dataPtr = dataPtr.getPtr(); }
  
  /// Checks validity of m_dataPtr.
  bool isMethodDataNull() const  {  return m_dataPtr.isNull(); }
  
  /// Gets the Class name
  static std::string getClassName() {  return DATA::getClassName() + "Strategy"; }

private: // data
  
  /// Pointer to the data object
  Common::SafePtr<DATA> m_dataPtr;
  
}; // end of class MethodStrategy

  /// This class represents a Null Strategy owned by a Method.
  /// @author Tiago Quintino
template <class DATA>
class Framework_API NullMethodStrategy : public MethodStrategy<DATA>    {
public:

  /// Constructor.
  explicit NullMethodStrategy(const std::string& name) : MethodStrategy<DATA>(name) {}

  /// Default destructor
  ~NullMethodStrategy() {}

  /// Checks if this object is a Null object.
  /// Since this is a Null command, it returns true
  /// @return true if Null and false otherwise
  bool isNull() const  {   return true;  }

}; // end of class NullMethodStrategy

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "MethodStrategy.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MethodStrategy_hh
