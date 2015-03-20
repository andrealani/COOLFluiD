// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshFormatConverter_hh
#define COOLFluiD_Framework_MeshFormatConverter_hh

//////////////////////////////////////////////////////////////////////////////


#include <boost/filesystem/path.hpp>

#include "Config/ConfigObject.hh"
#include "Common/OwnedObject.hh"
#include "Common/NonCopyable.hh"
#include "Common/NullableObject.hh"
#include "Common/FilesystemException.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ConcreteProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Provides an abstract interface for the format converters of the mesh files
/// in input.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API MeshFormatConverter : public Common::OwnedObject,
                            public Config::ConfigObject,
                            public Common::NonCopyable<MeshFormatConverter>,
                            public Common::NullableObject {
public: // types

  typedef Environment::ConcreteProvider<MeshFormatConverter,1> PROVIDER;
  typedef const std::string& ARG1;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  MeshFormatConverter(const std::string& name);

  /// Default destructor
  virtual ~MeshFormatConverter();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Tries to check the file for conformity to the format.
  /// Possibly not full proof.
  virtual void checkFormat(const boost::filesystem::path& filepath) = 0 ;

  /// Writes the data read to the original format.
  /// Useful for debugging purposes.
  virtual void convertBack(const boost::filesystem::path& filepath) = 0 ;

  /// Converts data from the file format to another format,
  /// taking into account the numerical method that will be
  /// used
  /// @param convertFromFileName name of the file to convert from
  /// @param fileName name of the file to write to
  virtual void convert(const boost::filesystem::path& fromFilepath,
                       const boost::filesystem::path& filepath);

  /// Gets the target file extention.
  virtual std::string getTargetExtension() const
  {
    return std::string("." + getTargetFormat());
  }

  /// Gets the origin file extension.
  virtual std::string getOriginExtension() const
  {
    return std::string("." + getOriginFormat());
  }
 
  /// Tell if this conveter can work in parallel.
  virtual bool isParallel() const {return false;}
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "MeshFormatConverter";
  }

protected: // functions

  /// Reads the all files needed for the conversion and assembles the data in the converter.
  /// It is always in a convert or convertBack, but only exectuted once.
  /// @throw Common::FilesystemException if a file cannot be open
  /// @throw BadFormatException if a file is ill formated
  virtual void readFiles(const boost::filesystem::path& filepath) = 0;

  /// Gets the flag with information if the solution
  /// space to be created should be discontinuous
  bool isDiscontinuous()
  {
    return m_isDiscontinuous;
  }

  /// Gets the target format.
  virtual std::string getTargetFormat() const = 0;

  /// Gets the origin format.
  virtual std::string getOriginFormat() const = 0;

  /// Gets the dimension.
  virtual CFuint getDimension() const
  {
    return Framework::PhysicalModelStack::getActive()->getDim();
  }

  /// Gets the number of variables
  virtual CFuint getNbVariables() const
  {
    return Framework::PhysicalModelStack::getActive()->getNbEq();
  }

  /// Adjust the node (state) numbering to make it stick to the
  /// COOLFluiD convention.
  virtual void adjustToCFmeshNodeNumbering() = 0;

  /// Write in the COOLFluiD format the element list for a FEM mesh
  virtual void writeContinuousElements(std::ofstream& fout) {}

  /// Write in the COOLFluiD format the state list for a FEM mesh
  virtual void writeContinuousStates(std::ofstream& fout) {}

  /// Write in the COOLFluiD format the Topological Region Set data
  /// considering to have FEM
  virtual void writeContinuousTrsData(std::ofstream& fout) {}

  /// Write in the COOLFluiD format the mesh data considering
  /// to have cell center FVM
  virtual void writeDiscontinuousElements(std::ofstream& fout) {}

  /// Write in the COOLFluiD format the Topological Region Set data
  /// considering to have FVM
  virtual void writeDiscontinuousTrsData(std::ofstream& fout) {}

  /// Write in the COOLFluiD format the state list for a
  /// cell centered FVM mesh
  virtual void writeDiscontinuousStates(std::ofstream& fout) {}

  /// Write in the COOLFluiD format the node list
  virtual void writeNodes(std::ofstream& fout) {}

private: // data

  /// is the solution space discontinuous,
  /// for example as in FiniteVolume or DiscontinuousGalerkin
  bool m_isDiscontinuous;

protected:

  /// list of boundary patches to ignore while converting
  std::vector<std::string> m_ignoreTRSNames;

}; // end MeshFormatConverter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(MeshFormatConverter)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshFormatConverter_hh
