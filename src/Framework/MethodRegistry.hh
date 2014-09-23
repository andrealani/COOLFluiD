// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MethodRegistry_hh
#define COOLFluiD_Framework_MethodRegistry_hh

//////////////////////////////////////////////////////////////////////////////

#include <set>

#include "Common/CFLog.hh"
#include "Common/SafePtr.hh"
#include "Common/NonCopyable.hh"
#include "Common/GeneralStorage.hh"
#include "Framework/QualifiedName.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Method;

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Registration facility for all the Methods
/// This class is a singleton object.
/// @author Tiago Quintino
class Framework_API MethodRegistry : public Common::NonCopyable<MethodRegistry> {

private: // typedefs

  typedef Common::GeneralStorage< Method > MtdStorage;

  typedef std::set<QualifiedName, QNameComp> QNameStorage;

  typedef Common::GeneralStorage< QNameStorage > MTypeStorage;

public: // methods for singleton

  /// @return the instance of this singleton
  static MethodRegistry& getInstance();

  /// Regists the Method on the Registry
  /// @param type is the abstract type of the method
  /// @param namesp is the namespace to which the of the method belongs
  /// @param name is the name of the method
  template < typename TYPE >
  void regist ( TYPE * ptr )
  {
    std::string type = TYPE::getClassName();
    insure_mstorage(type);
    QualifiedName qname  = ptr->QName();
    m_storebytype.getEntry(type)->insert(qname);
    m_storebyqname.addEntry(qname.str(), ptr);
  }

  /// Removes the Method from teh registry
  /// @param type is the abstract type of the method
  /// @param namesp is the namespace to which the of the method belongs
  /// @param name is the name of the method
  template < typename TYPE >
  void unregist ( TYPE * ptr )
  {
    std::string type = TYPE::getClassName();
    insure_mstorage(type);
    QualifiedName qname  = ptr->QName();
    m_storebyqname.removeEntry(qname.str());
    m_storebytype.getEntry(type)->erase(qname);
  }

  /// @return the Method that matches the name passed.
  /// @param name name of the Method
  /// @deprecated should not be used because it breaks if methods have same name across namespaces
  Common::SafePtr<Method> getMethod(const std::string& name);

  /// @return the Method that matches the name passed.
  /// @param qname QualifiedName of the Method
  Common::SafePtr<Method> getMethod(const QualifiedName& qname);

  /// @return the Method that matches the type, the namespace and the name passed.
  /// @param type is the abstract type of the method
  /// @param qname QualifiedName of the Method
  Common::SafePtr<Method> getMethod(const std::string& type, const QualifiedName& qname);

  /// @return the Method that matches the type, the namspace and the name passed.
  /// @param name QualifiedName of the Method
  template < typename TYPE >
  Common::SafePtr<TYPE> getMethod (const QualifiedName& qname)
  {
    std::string type = TYPE::getClassName();
    insure_mstorage(type);
    Method * ptr = m_storebyqname.getEntry(qname.str());
    Common::SafePtr<TYPE> rptr = dynamic_cast<TYPE*>(ptr);
    if (rptr.isNull())
      throw Common::FailedCastException (FromHere(),"Method " + qname.str() + "does not beong to type " + type);
    return rptr;
  }

  /// @return all the Method's that match the type and the namspace passed.
  /// @param type is a string describing the Method type
  /// @param namesp is the namespace to which the of the method belongs
  std::vector< Common::SafePtr<Method> >
  getAllMethods(const std::string& type, const std::string& namesp);

  /// @return all the Method's that match the type and the namspace passed.
  /// @param namesp is the namespace to which the of the method belongs
  template < typename TYPE >
  std::vector< Common::SafePtr<TYPE> > getAllMethods (const std::string& namesp = std::string())
  {
    std::string type = TYPE::getClassName();
    insure_mstorage(type);

    QNameStorage * qnst = m_storebytype.getEntry(type);
    std::vector< Common::SafePtr<TYPE> > all;
    typename QNameStorage::iterator itr = qnst->begin();
    for (; itr != qnst->end(); ++itr)
    {
      if ( namesp.empty() || (itr->getNamespace() == namesp) )
      {
        Method * ptr = m_storebyqname.getEntry(itr->str());
        Common::SafePtr<TYPE> rptr = dynamic_cast<TYPE*>(ptr);
        cf_assert(rptr.isNotNull());
        all.push_back(rptr);
        CFLog(DEBUG_MED, "Returning method " << rptr->getName() << "\n");
      }
    }
    return all;
  }

  /// @return all the Method's in the registry
  std::vector< Common::SafePtr<Method> > getAllMethods();

private: // helper functions

  /// Destructor removes all entries in storage
  ~MethodRegistry();

  std::string make_string_key(const std::string& type, const std::string& namesp);

  void insure_mstorage(const std::string& key);

private: // data

  /// storage of the QualifiedNames by Method TYPE
  MTypeStorage m_storebytype;
  /// storage of the pointe to the Method's associated to the QualifiedNames
  MtdStorage   m_storebyqname;

}; // end of class MethodRegistry

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Framework_MethodRegistry_hh
