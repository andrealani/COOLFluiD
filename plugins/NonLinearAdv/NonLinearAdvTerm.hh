// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_NonLinearAdv_NonLinearAdvTerm_hh
#define COOLFluiD_Physics_NonLinearAdv_NonLinearAdvTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NonLinearAdvTerm.
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class NonLinearAdvTerm : public Framework::BaseTerm {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {VX=0, VY=1, VZ=2};

  /**
   * Constructor without arguments
   */
  NonLinearAdvTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NonLinearAdvTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 3;
  }

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "NonLinearAdvTerm";
  }


private:


}; // end of class NonLinearAdvTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NonLinearAdv_NonLinearAdvTerm_hh
