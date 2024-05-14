// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CGNSWriter_CGWriterData_hh
#define COOLFluiD_IO_CGNSWriter_CGWriterData_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Framework/State.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/OutputFormatterData.hh"

#include "CGNSWriter/CGNSWriterAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class ConvectiveVarSet; class VarSetTransformer; }

    namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// CGNSWriterCom 's that compose the CGNSWriter.
/// @author Rayan Dhib

class CGNSWriter_API CGWriterData : public Framework::OutputFormatterData {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  CGWriterData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~CGWriterData();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Gets the Class name
  static std::string getClassName() { return "CGWriter"; }

  /// Sets the filename
  /// @param filepath pth to the file to be open
  void setFilename(const boost::filesystem::path& filepath)
  {
    m_filepath = filepath;
  }

  /// Gets the filename
  /// @return path to the file
  boost::filesystem::path getFilename() const
  {
    return m_filepath;
  }

  /// Informs if we should print extra values
  bool shouldPrintExtraValues() const
  {
    return m_printExtraValues;
  }

  /// outut only the surface file
  bool onlySurface() const
  {
    return m_surface_only;
  }

  /// Flag telling if to print only coordinates
  bool onlyCoordinates() const
  {
    return m_coord_only;
  }

  /// Flag telling if to print the equations
  bool withEquations() const
  {
    return m_withEquations;
  }
  
  /// Gets the filename
  /// @todo missing documentation
  Common::SafePtr<Framework::ConvectiveVarSet> getUpdateVarSet() const
  {
    return m_updateVarSet;
  }
  
  /// Gets the filename
  /// @todo missing documentation
  Common::SafePtr<Framework::ConvectiveVarSet> getOutputVarSet() const
  {
    return m_outputVarSet.getPtr();
  }
  
  /// Accessor to appendAuxData
  const bool getAppendAuxData() const
  {
    return m_appendAuxData;
  }

  Common::SafePtr<Framework::VarSetTransformer> getUpdateToOutputVarSetTransformer()
  {
    return m_updateToOutputVar.getPtr();
  }

  /// Gets the list of surface TRS names to write
  /// @return a vector of std::string with the names
  std::vector<std::string>& getSurfaceTRSsToWrite() {  return m_surfTRS; }

private:

  /// Filename to write solution to.
  boost::filesystem::path m_filepath;

  /// Name of the update variable set
  std::string m_outputVarStr;
  
  /// Flag telling if to print extra values
  bool m_printExtraValues;

  /// print only the surface data
  bool m_surface_only;

  /// print only the coordinates
  bool m_coord_only;

  /// print the equations
  bool m_withEquations;
  
  /// Update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_updateVarSet;

  /// Update variable set
  Common::SelfRegistPtr<Framework::ConvectiveVarSet> m_outputVarSet;

  /// list of names of TRS's to write in the surface file
  std::vector<std::string> m_surfTRS;
  
  /// Transformer from update to output Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_updateToOutputVar;

  /// Switch wether to write auxiliary data
  bool m_appendAuxData;

}; // end of class CGWriterData

//////////////////////////////////////////////////////////////////////////////

//  Comparision Funtion Object for two diferent state pointers
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

/// Definition of a command for CGWriter
typedef Framework::MethodCommand<CGWriterData> CGWriterCom;

/// Definition of a command provider for CGWriter
typedef Framework::MethodCommand<CGWriterData>::PROVIDER CGWriterComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CGNSWriter_CGWriterData_hh
