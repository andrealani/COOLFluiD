// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_MarcoTest_MarcoTestMethodData_hh
#define COOLFluiD_Numerics_MarcoTest_MarcoTestMethodData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvergenceMethodData.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceCellTrsGeoBuilder.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/GeoDataComputer.hh"
#include "Framework/NodalStatesExtrapolator.hh"
#include "MarcoTest/MarcoTestAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a Data Object that is accessed by the different
  /// MarcoTestCom 's that compose the MarcoTest.
  /// @see MarcoTestCom
  /// @author Andrea Lani
  /// @author Marco Panesi
class MarcoTest_API MarcoTestMethodData : public Framework::ConvergenceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  MarcoTestMethodData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~MarcoTestMethodData();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );
  
  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "MarcoTestMethod";
  }
  
  /// @return the GeometricEntity builder for faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
  getFaceTrsGeoBuilder()
  {
    return &m_faceTrsGeoBuilder;
  }

  /// @return the GeometricEntity builder that creates faces and cells on both sides (L/R)
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceCellTrsGeoBuilder> >
    getFaceCellTrsGeoBuilder()
  {
    return &m_faceCellTrsGeoBuilder;
  }
  
  /// @return the GeometricEntity builder for cell
  Common::SafePtr<Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> >
  getCellTrsGeoBuilder()
  {
    return &m_cellTrsGeoBuilder;
  }

  /// @return the GeometricEntity builder
  Common::SafePtr<Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> >
    getGeoWithNodesBuilder()
  {
    return &m_geoWithNodesBuilder;
  }
  
  /// Get the computer of geometric data
  /// @return a reference to the computer of geometric data
  Common::SafePtr<Framework::GeoDataComputer<MarcoTestMethodData> > getGeoDataComputer() const
  {
    cf_assert(m_geoDataComputer.isNotNull());
    return m_geoDataComputer.getPtr();
  }  
  
  /// Get the nodal extrapolator
  Common::SafePtr<Framework::NodalStatesExtrapolator<MarcoTestMethodData> >
    getNodalStatesExtrapolator() const
  {
    cf_assert(m_nStatesExtrapolator.isNotNull());
    return m_nStatesExtrapolator.getPtr();
  }
  
private:  
  
  /// Configures the GeoDataComputer
  void configureGeoDataComputer ( Config::ConfigArgs& args );
  
  /// Configures the NodalStatesExtrapolator
  void configureNodalStatesExtrapolator ( Config::ConfigArgs& args );
  
private: // data
  
  /// builder for faces
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceTrsGeoBuilder;
  
  /// builder for faces with cells
  Framework::GeometricEntityPool<Framework::FaceCellTrsGeoBuilder> m_faceCellTrsGeoBuilder;
  
  /// builder for cells
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_cellTrsGeoBuilder;
  
  // builder of GeometricEntity's with Node's
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> m_geoWithNodesBuilder;
  
  /// Geometric data computer
  Common::SelfRegistPtr<Framework::GeoDataComputer<MarcoTestMethodData> > m_geoDataComputer;
  
  /// nodal states extrapolator
  Common::SelfRegistPtr<Framework::NodalStatesExtrapolator<MarcoTestMethodData> > m_nStatesExtrapolator;
  
  /// string for the configuration of the geometric data computer
  std::string m_geoDataComputerStr;
  
  ///  name of the nodal extrapolator
  std::string m_nStatesExtrapolatorStr;
  
}; // end of class MarcoTestMethodData
    
//////////////////////////////////////////////////////////////////////////////
    
/// Definition of a command for MarcoTest
typedef Framework::MethodCommand<MarcoTestMethodData> MarcoTestMethodCom;

/// Definition of a command provider for MarcoTest
typedef Framework::MethodCommand<MarcoTestMethodData>::PROVIDER MarcoTestMethodComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace MarcoTest

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MarcoTest_MarcoTestMethodData_hh
