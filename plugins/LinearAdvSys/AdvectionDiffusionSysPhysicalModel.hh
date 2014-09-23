// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_AdvectionDiffusionSysPhysicalModel_hh
#define COOLFluiD_Physics_LinearAdvSys_AdvectionDiffusionSysPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "ADSysTerm.hh"
#include "LinearAdvSys/LinearAdvSysTerm.hh"
#include "Framework/ConvectionDiffusionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics {
namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a AdvectionDiffusionPhysicalModel.
/// @author Nadege Villedieu

template <int DIM>
class AdvectionDiffusionSysPhysicalModel :
public Framework::ConvectionDiffusionPM
<LinearAdvSysTerm, ADSysTerm> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

    /// Constructor without arguments
  AdvectionDiffusionSysPhysicalModel(const std::string& name);

  /// Default destructor
  virtual ~AdvectionDiffusionSysPhysicalModel();

  /// Configures this object by complementing the
  /// implementation in ConfigObject
  virtual void configure ( Config::ConfigArgs& args );

  /// Get the name of the Physical Model type
  /// @return name in a std::string
  std::string getTypeName() const
  {
    return std::string("AdvectionDiffusionSys" + Common::StringOps::to_str(DIM) + "D");
  }

  /// Get the convective name
  std::string getConvectiveName() const;

  /// Get the diffusive name
  std::string getDiffusiveName() const;

  /// @return the space dimension of the SubSystem
  virtual CFuint getDimension() const;

  /// @return the number of equations of the SubSystem
  virtual CFuint getNbEquations() const;


private:

  virtual void setReferenceValues()
    {
    }

  virtual void setReferenceTime()
  {
    /// @warning the reference time set below is only okay if the reference length is equal
    /// to one, or if the diffusion coefficient is equal to zero. To make it general, the diffusion
    /// coefficient should be rescaled also.
    _refTime = getRefLength();
  }


}; // end of class AdvectionDiffusionPhysicalModel

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

#include "AdvectionDiffusionSysPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_AdvectionDiffusionPhysicalModel_hh
