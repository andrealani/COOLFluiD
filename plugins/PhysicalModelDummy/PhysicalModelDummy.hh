// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_PhysicalModelDummy_PhysicalModelDummy_hh
#define COOLFluiD_PhysicalModelDummy_PhysicalModelDummy_hh

#include "Framework/ConvectionPM.hh"
#include "PhysicalModelDummy/DummyTerm.hh"

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a PhysicalModelDummy
class PhysicalModelDummy : public Framework::ConvectionPM< DummyTerm > {

public:

  /// Constructor without arguments
  PhysicalModelDummy(const std::string& name);

  /// Default destructor
  ~PhysicalModelDummy();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// @return name of the Physical Model type (std::string)
  std::string getTypeName() const
  {
    return std::string("PhysicalModelDummy");
  }

  /// @return the convective name
  std::string getConvectiveName() const
  {
    return getTypeName();
  }

  /// @return the diffusive name
  std::string getDiffusiveName() const
  {
    return "Null";
  }

  /// @return the space dimension of the SubSystem
  CFuint getDimension() const
  {
    return m_dim;
  }

  /// @return the number of equations of the SubSystem
  CFuint getNbEquations() const
  {
    return m_var;
  }

  /// Configures this object by complementing implementation in ConfigObject
  virtual void configure ( Config::ConfigArgs& args );


private:  // methods

  /// Set reference values
  void setReferenceValues();

  /**
   * Set the reference value for time
   */
  void setReferenceTime();

private:  // data

  /// Number of dimensions
  CFuint m_dim;

  /// Vector of equations' names
  std::vector< std::string > m_varnames;

  /// Number of equations per State
  CFuint m_var;

}; // end of class PhysicalModelDummy

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

#endif // COOLFluiD_PhysicalModelDummy_PhysicalModelDummy_hh

