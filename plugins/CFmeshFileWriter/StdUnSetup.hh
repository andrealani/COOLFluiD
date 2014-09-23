// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshFileWriter_StdUnSetup_hh
#define COOLFluiD_IO_CFmeshFileWriter_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "CFmeshFileWriter/CFmeshWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class CFmeshFileWriter_API StdUnSetup : public CFmeshWriterCom {
public:

  /// Constructor.
  explicit StdUnSetup(std::string name)
    : CFmeshWriterCom(name)
  {
  }

  /// Destructor.
  ~StdUnSetup()
  {
  }

  /// Execute Processing actions
  void execute();

}; // class StdUnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileWriter_StdUnSetup_hh

