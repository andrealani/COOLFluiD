// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshFileReader_StdUnSetup_hh
#define COOLFluiD_IO_CFmeshFileReader_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "CFmeshFileReader/CFmeshReaderData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to
  /// unsetup the data from the CFmeshFileReader method
class CFmeshFileReader_API StdUnSetup : public CFmeshReaderCom {
public: // functions

  /// Constructor.
  explicit StdUnSetup(const std::string& name);

  /// Virtual destructor.
  virtual ~StdUnSetup();

  /// Execute the unsetup action
  void execute();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileReader_StdUnSetup_hh

