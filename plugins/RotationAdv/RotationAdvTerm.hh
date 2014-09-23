// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdv_RotationAdvTerm_hh
#define COOLFluiD_Physics_RotationAdv_RotationAdvTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdv {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a RotationAdvTerm.
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class RotationAdvTerm : public Framework::BaseTerm {

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
  enum {u=0, VX=1, VY=2, VZ=3};

  /**
   * Constructor without arguments
   */
  RotationAdvTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RotationAdvTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 4;
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
    return "RotationAdvTerm";
  }

  /**
   * Get the origin in XX coordinates
   */
  CFreal getOX() const
  {
    return m_OX;
  }

  /**
   * Get the origin in YY coordinates
   */
  CFreal getOY() const
  {
    return m_OY;
  }

  /**
   * Get the origin in ZZ coordinates
   */
  CFreal getOZ() const
  {
    return m_OZ;
  }

  /**
   * Check if rotation is clockwise
   */
  bool isClockwise() const
  {
    return m_clockwise;
  }


private:

  /// Center of Rotation X coordinate
  CFreal m_OX;

  /// Center of Rotation Y coordinate
  CFreal m_OY;

  /// Center of Rotation Z coordinate
  CFreal m_OZ;

  /// Clockwise Rotation flag
  bool m_clockwise;

}; // end of class RotationAdvTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_RotationAdv_RotationAdvTerm_hh
