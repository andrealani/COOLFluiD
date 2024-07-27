// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CGNSWriter_StdSetup_hh
#define COOLFluiD_IO_CGNSWriter_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "CGNSWriter/CGWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class CGNSWriter_API StdSetup : public CGWriterCom {
public:

  /// Constructor.
  explicit StdSetup(std::string name);

  /// Destructor.
  ~StdSetup();

  /// Execute Processing actions
  void execute();

}; // class StdSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CGNSWriter_StdSetup_hh

