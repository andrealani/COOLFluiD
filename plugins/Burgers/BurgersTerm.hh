// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_Burgers_BurgersTerm_hh
#define COOLFluiD_Physics_Burgers_BurgersTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Burgers {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a BurgersTerm.
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */
class BurgersTerm : public Framework::BaseTerm {

public:

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum {VX=0, VY=1};

  /**
   * Constructor without arguments
   */
  BurgersTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~BurgersTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 2;
  }

  /**
   * Get the name
   */
  static std::string getName()
  {
    return "BurgersTerm";
  }

}; // end of class BurgersTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Burgers

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Burgers_BurgersTerm_hh
