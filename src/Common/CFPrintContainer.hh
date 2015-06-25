// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CFPrintContainer_hh
#define COOLFluiD_Common_CFPrintContainer_hh

//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "Common/COOLFluiD.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class prints the entries in a given container via its overloaded operator <<
/// @author Andrea Lani
template <typename ARRAY>
class CFPrintContainer {

public:

  /// Constructor
  /// @param sequence name
  CFPrintContainer(std::string name, ARRAY *const v, CFuint stride = 1) :
    m_name(name), m_v(v), m_stride(stride)
  {
    cf_assert(m_stride > 0);
  }

  /// Constructor
  /// @param sequence name
  CFPrintContainer(const char name[], ARRAY *const v, CFuint stride = 1) :
    m_name(name), m_v(v), m_stride(stride)
  {
  }

  /// Default destructor
  ~CFPrintContainer()
  {
  }

  /// Overloading of streaming operator
  friend std::ostream& operator<< (std::ostream& out, const CFPrintContainer<ARRAY>& a)
  {
    out << "P" << Common::PE::GetPE().GetRank("Default") << " => " << a.m_name << " ";
    CFuint n = 0;
    typename ARRAY::const_iterator it;
    for (it = a.m_v->begin(); it != a.m_v->end(); ++it, ++n) {
      if ((n+1)%a.m_stride > 0) {
	out << *it << " ";
      }
      else {
	out << *it << "  ";
      }
    }
    out << "\n";
    return out;
  }
  
private:

  // message to print
  std::string m_name;

  // container to print
  ARRAY *const m_v;

  // stride in printing
  CFuint m_stride;

}; // end of class CFPrintContainer

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CFPrintContainer_hh
