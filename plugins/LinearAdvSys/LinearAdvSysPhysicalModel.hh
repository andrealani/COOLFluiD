// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_LinearAdvSysPhysicalModel_hh
#define COOLFluiD_Physics_LinearAdvSys_LinearAdvSysPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "LinearAdvSysTerm.hh"
#include "Framework/ConvectionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a LinearAdvSysPhysicalModel.
 *
 * @author Tiago Quintino
 * @author Tomas Kopacek
 *
 */

template <int DIM>
class LinearAdvSysPhysicalModel : public Framework::ConvectionPM<LinearAdvSysTerm> {
public:

  /**
   * Constructor without arguments
   */
  LinearAdvSysPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~LinearAdvSysPhysicalModel();

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("LinearAdvSys" + Common::StringOps::to_str(DIM) + "D");
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

}; // end of class LinearAdvSysPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

#include "LinearAdvSysPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdvSys_LinearAdvSysPhysicalModel_hh
