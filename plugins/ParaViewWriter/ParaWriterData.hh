// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_ParaViewWriter_ParaWriterData_hh
#define COOLFluiD_IO_ParaViewWriter_ParaWriterData_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Framework/GeometricEntityPool.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/OutputFormatterData.hh"
#include "Framework/State.hh"
#include "Framework/StdTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class ConvectiveVarSet;
  }

  namespace IO {

    namespace ParaViewWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Data Object that is accessed by the different
 * ParaViewWriterCom 's that compose the ParaViewWriter.
 *
 * @see ParaWriterCom
 *
 * @author Kris Van den Abeele
 *
 */
class ParaWriterData : public Framework::OutputFormatterData {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  ParaWriterData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~ParaWriterData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ParaWriter";
  }

  /**
   * Sets the filename
   * @param filepath pth to the file to be open
   */
  void setFilename(const boost::filesystem::path& filepath)
  {
    m_filepath = filepath;
  }

  /**
   * Gets the filename
   * @return path to the file
   */
  boost::filesystem::path getFilename() const
  {
    return m_filepath;
  }

  /**
   * Tells if to print extra values
   */
  bool printExtraValues() const
  {
    return m_printExtraValues;
  }

  /// ouptut only the surface file
  bool onlySurface() const
  {
    return m_surface_only;
  }

  /**
   * Gets the UpdateVarSet
   * @todo missing documentation
   */
  Common::SafePtr<Framework::ConvectiveVarSet> getUpdateVarSet() const
  {
    return m_updateVarSet.getPtr();
  }

  /**
   * Gets the list of surface TRS names to write
   * @return a vector of std::string with the names
   */
  std::vector<std::string>& getSurfaceTRSsToWrite()
  {
    return m_surfTRS;
  }

  /// @return the GeometricEntity builder
  Common::SafePtr<
      Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder > >
      getStdTrsGeoBuilder()
  {
    return &m_stdTrsGeoBuilder;
  }

  /**
   * returns the VTK cell type ID
   */
  static CFuint getVTKCellTypeID(CFGeoShape::Type shape,CFuint geoOrder);

  /**
   * returns the mapped coordinates of the output points in a cell of given shape and order
   */
  static std::vector< RealVector > getOutputPntsMappedCoords(CFGeoShape::Type shape,CFuint solOrder);

  /**
   * returns the local subcell node connectivity in a cell of given shape and order
   */
  static std::vector< std::vector< CFuint > > getOutputCellNodeConn(CFGeoShape::Type shape,CFuint solOrder);

  /// Accessor to the switch to write velocity by components or coupled
  bool writeVectorAsComponents() const
  {
    return m_writeVectorAsComponents;
  }

private:

  /// Filename to write solution to.
  boost::filesystem::path m_filepath;

  /// Name of the update variable set
  std::string m_updateVarStr;

  /// Flag telling if to print extra values
  bool m_printExtraValues;

  /// print only the surface data
  bool m_surface_only;

  /// Update variable set
  Common::SelfRegistPtr<Framework::ConvectiveVarSet> m_updateVarSet;

  /// list of names of TRS's to write in the surface file
  std::vector<std::string> m_surfTRS;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder >  m_stdTrsGeoBuilder;

  /// Switch to write velocity by components or coupled
  bool m_writeVectorAsComponents;

}; // end of class ParaWriterData

//////////////////////////////////////////////////////////////////////////////

/**
  Comparision Funtion Object for two diferent state pointers
 */
template <class TYPE>
struct LessLocalIDX : public std::binary_function<TYPE*,
                                                  TYPE*,
                                                  bool> {
  bool operator()(TYPE* x, TYPE* y)
  {
    return x->getLocalID() < y->getLocalID();
  }

 }; // end struct LessLocalIDX

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for ParaWriter
typedef Framework::MethodCommand<ParaWriterData> ParaWriterCom;

/// Definition of a command provider for ParaWriter
typedef Framework::MethodCommand<ParaWriterData>::PROVIDER ParaWriterComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_ParaViewWriter_ParaWriterData_hh
