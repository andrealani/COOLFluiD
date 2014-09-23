// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdvSys_RotationAdvSysPhysicalModel_hh
#define COOLFluiD_Physics_RotationAdvSys_RotationAdvSysPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionPM.hh"
#include "RotationAdvSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a RotationAdvPhysicalModel.
 *
 * @author Tiago Quintino
 *
 */

template <int DIM>
class RotationAdvSysPhysicalModel : public Framework::ConvectionPM<RotationAdvSysTerm> {

public:

  /**
   * Constructor without arguments
   */
  RotationAdvSysPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~RotationAdvSysPhysicalModel();

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("RotationAdvSys" + Common::StringOps::to_str(DIM) + "D");
  }

  /**
   * Get the convective name
   */
  std::string getConvectiveName() const
  {
    return getTypeName();
  }

  /**
   * Get the diffusive name
   */
  std::string getDiffusiveName() const
  {
    return "Null";
  }

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

private:

  /**
   * Set the reference values
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceValues();
  

  /**
   * Set the reference value for time
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceTime();

}; // end of class RotationAdvPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

#include "RotationAdvSysPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_RotationAdv_RotationAdvPhysicalModel_hh
