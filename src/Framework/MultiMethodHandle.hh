// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MultiMethodHandle_hh
#define COOLFluiD_Framework_MultiMethodHandle_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class holds multi method-related storage, self-registration
/// keys and names, allowing to perform configuration actions uniformly
/// for the different arrays of methods
/// @author Andrea Lani
template <class TMETHOD>
class  MultiMethodHandle {
public:

  /// Constructor
  MultiMethodHandle() : _mList() { }

  /// Destructor
  ~MultiMethodHandle() { }

  /// Copy Constructor
  MultiMethodHandle(const MultiMethodHandle<TMETHOD>& other) :  _mList(other._mList)  {}

  /// Assignment operator
  const MultiMethodHandle& operator=
    (const MultiMethodHandle<TMETHOD>& other)
  {
    _mList = other._mList;
    return *this;
  }

  /// Clear the object
  void clear() { std::vector<TMETHOD*>().swap(_mList); }
  
  /// Add method in the list
  void addMethod(TMETHOD *const method) { _mList.push_back(method);  }
  
  /// Apply given void action (pointer to member function with 0
  /// arguments) to all the methods in the list.
  void apply (std::mem_fun_t<void,TMETHOD> fun, CFint ns=-1)
  {
    const int rank = Common::PE::GetPE().GetRank("Default");
    const CFint nsize =  (ns < 0) ? _mList.size() : ns;
    for (CFuint i = 0; i < nsize; ++i) {
      if (Common::PE::GetPE().isRankInGroup(rank, _mList[i]->getNamespace())) {
        fun(_mList[i]);
      }
    }
  }
  
  /// Apply given void action (pointer to member function with 1
  /// argument) to all the methods in the list.
  template <class ARG>
  void apply (std::mem_fun1_t<void, TMETHOD, ARG> fun, ARG arg, CFint ns=-1)
  {
    const int rank = Common::PE::GetPE().GetRank("Default");
    const CFint nsize =  (ns < 0) ? _mList.size() : ns;
    for (CFuint i = 0; i < nsize; ++i) {
      if (Common::PE::GetPE().isRankInGroup(rank, _mList[i]->getNamespace())) {
	(std::bind2nd<std::mem_fun1_t<void, TMETHOD, ARG> >
	 (fun, arg))(_mList[i]);
      }
    }
  }
  
  /// Overloading of the idx operator
  Common::SafePtr<TMETHOD> operator[] (CFuint i) const
  {
    cf_assert(i < _mList.size());
    return _mList[i];
  }

  /// Get the size
  CFuint size() const { return _mList.size(); }

  /// Tell if the handle is null
  bool isNotNull() const { return (_mList.size() > 0); }

private:

  /// method list for a quick indexable access
  std::vector<TMETHOD*> _mList;

}; // end class MultiMethodHandle

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MultiMethodHandle_hh
