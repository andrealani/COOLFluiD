// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshFileWriter_StdSetup_hh
#define COOLFluiD_IO_CFmeshFileWriter_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "CFmeshFileWriter/CFmeshWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class CFmeshFileWriter_API StdSetup : public CFmeshWriterCom {
public:

  /// Constructor.
  explicit StdSetup(std::string name)
    : CFmeshWriterCom(name)
  {
  }

  /// Destructor.
  ~StdSetup()
  {
  }

  /// Execute Processing actions
  void execute();

}; // class StdSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileWriter_StdSetup_hh

