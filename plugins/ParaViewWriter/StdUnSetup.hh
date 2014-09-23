// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_ParaViewWriter_StdUnSetup_hh
#define COOLFluiD_IO_ParaViewWriter_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "ParaWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace ParaViewWriter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   */
class StdUnSetup : public ParaWriterCom {
public:

  /**
   * Constructor.
   */
  explicit StdUnSetup(std::string name)
    : ParaWriterCom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdUnSetup()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

}; // class StdUnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_ParaViewWriter_StdUnSetup_hh

