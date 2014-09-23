// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdvSys_RotationAdvSysTerm_hh
#define COOLFluiD_Physics_RotationAdvSys_RotationAdvSysTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a RotationAdvTerm.
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class RotationAdvSysTerm : public Framework::BaseTerm {

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
    enum {C0X=0, C0Y=1, C0Z=2, C1X=3, C1Y=4, C1Z=5, C2X=6, C2Y=7, C2Z=8, C3X=9, C3Y=10, C3Z=11, u0=12, u1=13, u2=14, u3=15};
  /**
   * Constructor without arguments
   */
  RotationAdvSysTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RotationAdvSysTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

/**
 *    * Resize the physical data
 *       */
  virtual void resizePhysicalData(RealVector& physicalData);

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 16;
  }

  /**
 *    * Get the start of the scalar vars data (mass fractions)
 *       */
  CFuint getFirstScalarVar(CFuint i) const
  {
    return 12;
  }

  /**
 *    * Get the number of scalar vars
 *       */
  CFuint getNbScalarVars(CFuint i) const
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
    return "RotationAdvSysTerm";
  }

  /**
   * Get the origin in XX coordinates
   */
  CFreal getOX0() const
  {
    return m_OX0;
  }

  /**
   * Get the origin in YY coordinates
   */
  CFreal getOY0() const
  {
    return m_OY0;
  }

  /**
   * Get the origin in ZZ coordinates
   */
  CFreal getOZ0() const
  {
    return m_OZ0;
  }


  /**
 *    * Get the origin in XX coordinates
 *       */
  CFreal getOX1() const
  {
    return m_OX1;
  }

  /**
 *    * Get the origin in YY coordinates
 *       */
  CFreal getOY1() const
  {
    return m_OY1;
  }

  /**
 *    * Get the origin in ZZ coordinates
 *       */
  CFreal getOZ1() const
  {
    return m_OZ1;
  }


  /**
 *    * Get the origin in XX coordinates
 *       */
  CFreal getOX2() const
  {
    return m_OX2;
  }

  /**
 *    * Get the origin in YY coordinates
 *       */
  CFreal getOY2() const
  {
    return m_OY2;
  }

  /**
 *    * Get the origin in ZZ coordinates
 *       */
  CFreal getOZ2() const
  {
    return m_OZ2;
  }

  /**
 *    * Get the origin in XX coordinates
 *       */
  CFreal getOX3() const
  {
    return m_OX3;
  }

  /**
 *    * Get the origin in YY coordinates
 *       */
  CFreal getOY3() const
  {
    return m_OY3;
  }

  /**
 *    * Get the origin in ZZ coordinates
 *       */
  CFreal getOZ3() const
  {
    return m_OZ3;
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
  CFreal m_OX0;

  /// Center of Rotation Y coordinate
  CFreal m_OY0;

  /// Center of Rotation Z coordinate
  CFreal m_OZ0;

  /// Center of Rotation X coordinate
  CFreal m_OX1;
  
  /// Center of Rotation Y coordinate
  CFreal m_OY1;
  
  /// Center of Rotation Z coordinate
  CFreal m_OZ1;
  
  
  /// Center of Rotation X coordinate
  CFreal m_OX2;
  
  /// Center of Rotation Y coordinate
  CFreal m_OY2;

    /// Center of Rotation Z coordinate
  CFreal m_OZ2;

  /// Center of Rotation X coordinate
  CFreal m_OX3;
      
  /// Center of Rotation Y coordinate
  CFreal m_OY3;
  
  /// Center of Rotation Z coordinate
  CFreal m_OZ3;
  
  
  /// Clockwise Rotation flag
  bool m_clockwise;

}; // end of class RotationAdvTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace RotationAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_RotationAdv_RotationAdvTerm_hh
