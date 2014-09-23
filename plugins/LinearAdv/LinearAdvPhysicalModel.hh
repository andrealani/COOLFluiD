// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdv_LinearAdvPhysicalModel_hh
#define COOLFluiD_Physics_LinearAdv_LinearAdvPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////
#include "LinearAdvTerm.hh"
#include "Framework/ConvectionPM.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a LinearAdvPhysicalModel.
/// @author Tiago Quintino
template <int DIM>
class LinearAdvPhysicalModel : public Framework::ConvectionPM<LinearAdvTerm> {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor without arguments
  LinearAdvPhysicalModel(const std::string& name);

  /// Default destructor
  ~LinearAdvPhysicalModel();

  /// Get the name of the Physical Model type
  /// @return name in a std::string
  std::string getTypeName() const
  {
    return std::string("LinearAdv" + Common::StringOps::to_str(DIM) + "D");
  }

  /// Get the convective name
  std::string getConvectiveName() const
  {
    return getTypeName();
  }

  /// Get the diffusive name
  std::string getDiffusiveName() const
  {
    return "Null";
  }

  /// @return the space dimension of the SubSystem
  CFuint getDimension() const;

  /// @return the number of equations of the SubSystem
  CFuint getNbEquations() const;

  /// Configures this object by complementing the
  /// implementation in ConfigObject
  virtual void configure ( Config::ConfigArgs& args );

private:

  /// Set the reference values
  /// This function is virtual to allow to use a different
  /// adimensionalizaton, if needed
  virtual void setReferenceValues()
  {
  }

  /// Set the reference value for time
  /// This function is virtual to allow to use a different
  /// adimensionalizaton, if needed
  virtual void setReferenceTime()
  {
    /// @warning the reference time set below is only okay if the reference length is equal
    /// to one, or if the diffusion coefficient is equal to zero. To make it general, the diffusion
    /// coefficient should be rescaled also.
    _refTime = getRefLength();
  }

private:

  /// Advection speed in X direction
  CFreal m_VX;

  /// Advection speed in Y direction
  CFreal m_VY;

  /// Advection speed in Y direction
  CFreal m_VZ;

}; // end of class LinearAdvPhysicalModel

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "LinearAdvPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_LinearAdvPhysicalModel_hh
