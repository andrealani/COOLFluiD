// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataProcessingData_hh
#define COOLFluiD_Framework_DataProcessingData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/MethodData.hh"
#include "Framework/SpaceMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class ConvectiveVarSet;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// DataProcessingCom 's that compose the DataProcessingMethod
/// @todo there is missing documentation in this class.
/// @author Thomas Wuilbaut
class Framework_API DataProcessingData : public Framework::MethodData {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments.
  DataProcessingData(Common::SafePtr<Method> owner);

  /// Destructor.
  ~DataProcessingData();

  /// Configure the data from the supplied arguments.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Gets the Class name
  static std::string getClassName() {  return "DataProcessing"; }

  /// Gets the update var set
  Common::SafePtr<Framework::ConvectiveVarSet> getUpdateVarSet() const
  {
    return _updateVarSet.getPtr();
  }
  
  std::string getUpdateVarSetStr() const
  {
    return _updateVarStr;
  }

  /// @return the GeometricEntity builder
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> >
  getStdTrsGeoBuilder()
  {
    return &_stdTrsGeoBuilder;
  }

  /// @return the GeometricEntity builder for faces in FVMCC
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
  getFaceTrsGeoBuilder()
  {
    return &_faceTrsGeoBuilder;
  }
  
  /// get the iteration number at which processing should start
  CFuint getStartIter() const {return m_startIter;}
  
private:

  /// handle to the space method
  Framework::MultiMethodHandle<Framework::SpaceMethod> _spaceMtd;

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;

  /// builder for face TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceTrsGeoBuilder;

  /// Name of the update variable set
  std::string _updateVarStr;

  /// Update variable set
  Common::SelfRegistPtr<Framework::ConvectiveVarSet> _updateVarSet;

  /// iteration at which postprocessing starts
  CFuint m_startIter;
  
}; // end of class DataProcessingData

//////////////////////////////////////////////////////////////////////////////

  /// Definition of a command for DataProcessingMethod
  typedef Framework::MethodCommand<DataProcessingData> DataProcessingCom;

  /// Definition of a command provider for DataProcessingMethod
  typedef Framework::MethodCommand<DataProcessingData>::PROVIDER DataProcessingComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace DataProcessing

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataProcessingData_hh
