// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_TecplotWriter_WriteSolution1D_hh
#define COOLFluiD_TecplotWriter_WriteSolution1D_hh

//////////////////////////////////////////////////////////////////////////////

#include "TecplotWriter/WriteSolution.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to write the MeshData
/// solution in a 1D case to a Tecplot format  file for visualization.
/// @author Andrea Lani
class TecplotWriter_API WriteSolution1D : public WriteSolution {
public:

  /// Constructor.
  explicit WriteSolution1D(const std::string& name);

  /// Destructor.
  ~WriteSolution1D()
  {
  }

  /// Set up private data
  virtual void setup();

  /// Execute Processing actions
  void execute();

protected:

  /// Write the to the given file stream the MeshData.
  /// @throw Common::FilesystemException
  virtual void writeToFileStream(std::ofstream& fout);

}; // class WriteSolution1D

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_TecplotWriter_WriteSolution1D_hh

