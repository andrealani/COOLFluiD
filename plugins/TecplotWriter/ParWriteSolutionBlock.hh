// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_TecplotWriter_ParWriteSolutionBlock_hh
#define COOLFluiD_IO_TecplotWriter_ParWriteSolutionBlock_hh

//////////////////////////////////////////////////////////////////////////////

#include "TecplotWriter/ParWriteSolution.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a parallel writer for TECPLOT files in block format
/// @author Andrea Lani
class TecplotWriter_API ParWriteSolutionBlock : public ParWriteSolution {
  
public:
  
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Constructor.
  explicit ParWriteSolutionBlock(const std::string& name);

  /// Destructor.
  virtual ~ParWriteSolutionBlock();
  
  /// Set up private data
  virtual void setup();
  
 protected:
    
  /// Writes the boundary data to file
  /// @param filename  name of output file
  /// @param isNewFile flag to tell if the file is to be created or to be overwritten
  /// @param fout      pointer to the file  
  virtual void writeBoundaryData(const boost::filesystem::path& filepath,
				 const bool isNewFile,
				 const std::string title,
				 std::ofstream* fout);
  
  /// Writes the TECPLOT zone header
  virtual void writeZoneHeader(std::ofstream* fout, 
			       const CFuint iType,
			       const std::string& geoShape,
			       const CFuint nbNodesInType,
			       const CFuint nbElemsInType,
			       const std::string& geoType,
			       const std::string& end,
			       const bool isBoundary); 
  
  /// Write the Tecplot file in Binmary format
  /// @throw Common::FilesystemException
  virtual void writeToBinaryFile();
  
  /// Write the node list corresponding to the given element type
  virtual void writeNodeList(std::ofstream* fout, const CFuint iType, 
			     Common::SafePtr<Framework::TopologicalRegionSet> elements,
			     const bool isBoundary);
  
 protected:
  
  /// flag that specifies to output cell-centered or nodal variables
  bool m_nodalOutputVar;
  
}; // class ParWriteSolutionBlock

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_TecplotWriter_ParWriteSolutionBlock_hh

