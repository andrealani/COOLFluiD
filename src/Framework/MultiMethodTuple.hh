// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MultiMethodTuple_hh
#define COOLFluiD_Framework_MultiMethodTuple_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/CFMap.hh"
#include "Common/SelfRegistPtr.hh"
#include "Framework/Method.hh"
#include "Framework/MethodRegistry.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class holds multi method-related storage, self-registration
/// keys and names, allowing to perform configuration actions uniformly
/// for the different arrays of methods
/// @author Andrea Lani
template <class TMETHOD>
class MultiMethodTuple {
public: // functions

  /// Constructor
  MultiMethodTuple() : mKeys(), mNames(), m_mList(), m_mMap()
  {
  }

  /// Destructor
  ~MultiMethodTuple()
  {
    typename std::vector<TMETHOD*>::iterator itr = m_mList.begin();
    for (; itr != m_mList.end(); ++itr)
    {
      TMETHOD* ptr = *itr;
      try
      {
        MethodRegistry::getInstance().template unregist<TMETHOD>( ptr );
      }
      catch (Common::NoSuchStorageException& e)
      {
  CFLog(VERBOSE, e.what());
        CFLog(WARN, "Could not unregist method " + ptr->getName() + " of type " + TMETHOD::getClassName());
      }
    }
  }

  /// Apply given void action (pointer to member function with 0
  /// arguments) to the methods in the list that satify the given condition
  void apply (std::mem_fun_t<void,TMETHOD> fun, 
	      std::mem_fun_t<bool,TMETHOD> pr, 
	      bool flag, CFint ns=-1)
  {
    const int rank = Common::PE::GetPE().GetRank("Default");
    const CFint nsize =  (ns < 0) ? m_mList.size() : ns;
    for (CFuint i = 0; i < nsize; ++i) {
      if (Common::PE::GetPE().isRankInGroup(rank, m_mList[i]->getNamespace())) {
	if (pr(m_mList[i]) == flag) {
	  CFLog(VERBOSE, "MultiMethodTuple::apply() => rank [" << rank << 
		"], group ["<< m_mList[i]->getNamespace() << "], method [" << m_mList[i]->getName() << "]\n");
	  fun(m_mList[i]);
	}
      }
    }
  }
  
  /// Apply given void action (pointer to member function with 0
  /// arguments) to all the methods in the list.
  void apply (std::mem_fun_t<void,TMETHOD> fun, CFint ns=-1)
  {
    const int rank = Common::PE::GetPE().GetRank("Default");
    const CFint nsize =  (ns < 0) ? m_mList.size() : ns;
    for (CFuint i = 0; i < nsize; ++i) {
      if (Common::PE::GetPE().isRankInGroup(rank, m_mList[i]->getNamespace())) {
         CFLog(VERBOSE, "MultiMethodTuple::apply() => rank [" << rank << 
                "], group ["<< m_mList[i]->getNamespace() << "], method [" << m_mList[i]->getName() << "]\n");
        fun(m_mList[i]);
      }
    }
  }
  
  /// Apply given void action (pointer to member function with 1
  /// argument) to all the methods in the list.
  template <class ARG>
  void apply (std::mem_fun1_t<void, TMETHOD, ARG> fun, ARG arg, CFint ns=-1)
  {
    const int rank = Common::PE::GetPE().GetRank("Default");
    const CFint nsize =  (ns < 0) ? m_mList.size() : ns;
    for (CFuint i = 0; i < nsize; ++i) {
      if (Common::PE::GetPE().isRankInGroup(rank, m_mList[i]->getNamespace())) {
         CFLog(VERBOSE, "MultiMethodTuple::apply() => rank [" << rank << 
                "], group ["<< m_mList[i]->getNamespace() << "], method [" << m_mList[i]->getName() << "]\n");
        (std::bind2nd<std::mem_fun1_t<void, TMETHOD, ARG> >(fun, arg)) (m_mList[i]);
      }
    }  
  }
  
  /// Apply given void action (pointer to member function with 1
  /// argument) to the methods in the list that satify the given condition
  template <class ARG>
  void apply (std::mem_fun1_t<void, TMETHOD, ARG> fun, ARG arg, 
	      std::mem_fun_t<bool,TMETHOD> pr, 
	      bool flag, CFint ns=-1)
  {
    const int rank = Common::PE::GetPE().GetRank("Default");
    const CFint nsize =  (ns < 0) ? m_mList.size() : ns;
    for (CFuint i = 0; i < nsize; ++i) {
      if (Common::PE::GetPE().isRankInGroup(rank, m_mList[i]->getNamespace())) {
	if ( pr(m_mList[i]) == flag) {
	  CFLog(VERBOSE, "MultiMethodTuple::apply() => rank [" << rank << 
		"], group ["<< m_mList[i]->getNamespace() << "], method [" << m_mList[i]->getName() << "]\n");
	  (std::bind2nd<std::mem_fun1_t<void, TMETHOD, ARG> >(fun, arg)) (m_mList[i]);
	}
      }
    }  
  }
  
  /// Apply given void action (pointer to member function with 0
  /// arguments) to all the methods in the list.
  void apply (Method::root_mem_fun_t<void,TMETHOD> fun, CFint ns=-1)
  {
    const int rank = Common::PE::GetPE().GetRank("Default");
    const CFint nsize =  (ns < 0) ? m_mList.size() : ns;
    for (CFuint i = 0; i < nsize; ++i) {
      if (Common::PE::GetPE().isRankInGroup(rank, m_mList[i]->getNamespace())) {
         CFLog(VERBOSE, "MultiMethodTuple::apply() => rank [" << rank << 
                "], group ["<< m_mList[i]->getNamespace() << "], method [" << m_mList[i]->getName() << "]\n");
        fun(m_mList[i]);
      }
    }
  }
  
  /// Overloading of the idx operator
  TMETHOD* operator[] (CFuint i) const
  {
    return m_mList[i];
  }

  /// Find name and return corresponding method
  /// @throw Common::NoSuchValueException if the Method is not present in the Tuple
  TMETHOD* find (const std::string& name)
  {
    return m_mMap.find(name).getPtr();
  }


  /// Get the size
  CFuint size() const
  {
    cf_assert(m_mList.size() == m_mMap.size());
    return m_mList.size();
  }
  
  /// Get the beginning
  TMETHOD** begin()
  {
    cf_assert(m_mList.size() == m_mMap.size());
    cf_assert(m_mList.size() > 0);
    return &m_mList[0];
  }
  
  /// Get the end
  TMETHOD** end()
  {
    cf_assert(m_mList.size() == m_mMap.size());
    cf_assert(m_mList.size() > 0);
    return &m_mList[0] + m_mList.size();
  }
  
  /// Add a pointer to a new method to the list
  void addPtr(const std::string& name, Common::SelfRegistPtr<TMETHOD>& met)
  {
    m_mMap.insert(name,met);
    m_mList.push_back(met.getPtr());
  }

public: // data

  ///  Self-registering keys of Methods to instantiate
  std::vector<std::string> mKeys;

  /// Names of Methods to configure
  std::vector<std::string> mNames;

private:

  /// method list for a quick indexable access
  std::vector<TMETHOD*> m_mList;

  /// method map
  Common::CFMap<std::string, Common::SelfRegistPtr<TMETHOD> > m_mMap;

}; // end class MultiMethodTuple

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MultiMethodTuple_hh
