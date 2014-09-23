// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_NonLinearAdv_NonLinearAdvPhysicalModel_hh
#define COOLFluiD_Physics_NonLinearAdv_NonLinearAdvPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionPM.hh"
#include "NonLinearAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NonLinearAdvPhysicalModel.
 *
 * @author Nadege Villedieu
 *
 */

template <int DIM>
class NonLinearAdvPhysicalModel : public Framework::ConvectionPM<NonLinearAdvTerm> {

public:

  /**
   * Constructor without arguments
   */
  NonLinearAdvPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~NonLinearAdvPhysicalModel();

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("NonLinearAdv" + Common::StringOps::to_str(DIM) + "D");
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
  virtual void setReferenceValues()
  {
  }
  /**
   * Set the reference value for time
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceTime()
{
}

}; // end of class NonLinearAdvPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

#include "NonLinearAdvPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NonLinearAdv_NonLinearAdvPhysicalModel_hh
