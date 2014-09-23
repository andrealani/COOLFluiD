// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_CGNSException_hh
#define COOLFluiD_CGNS_CGNSException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "CGNS2CFmesh/CGNSDefinitions.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/**
 * Excepion thrown when problems happen when reading a CGNS file
 *
 * @author Tiago Quintino
 */
class CGNSException : public Common::Exception {
public:

  /**
   * Constructor
   *
   * @param what is the name of the file that has been requested, but
   *                  cannot actually be opened
   * @see Exception()
   */
  CGNSException (const Common::CodeLocation& where, const std::string& what) : Exception(where,what,"CGNSException")
  {
  }

  /**
   * A copy constructor is necessary for exceptions, for the C++
   * exception mechanism to work.
   */
  CGNSException(const CGNSException& e) throw() : Exception(e)
  {
  }

}; // end of class CGNSException

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_CGNSException_hh
