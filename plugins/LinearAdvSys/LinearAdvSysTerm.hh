// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_LinearAdvSysTerm_hh
#define COOLFluiD_Physics_LinearAdvSys_LinearAdvSysTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an Linear Advection for Systems
 * physical term
 *
 * @author Tiago Quintino
 * @author Tomas Kopacek
 *
 */
class LinearAdvSysTerm : public Framework::BaseTerm {
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
  LinearAdvSysTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LinearAdvSysTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Resize the physical data
   */
  virtual void resizePhysicalData(RealVector& physicalData);

  /**
   * Physical data size 
   */
  virtual CFuint getDataSize() const 
  {  
    return 16; 
  }

  /**
   * Get the start of the scalar vars data (mass fractions)
   */
  CFuint getFirstScalarVar(CFuint i) const
  {
    return 12;
  }

  /**
   * Get the number of scalar vars
   */
  CFuint getNbScalarVars(CFuint i) const
  {
    return 0;
  }
  
  /**
   * Get the convecting velocities
   */
  CFreal getc0x() const
  {
    return m_c0x;
  }

  CFreal getc0y() const
  {
    return m_c0y;
  }

  CFreal getc0z() const
  {
    return m_c0z;
  }

  CFreal getc1x() const
  {
    return m_c1x;
  }

  CFreal getc1y() const
  {
    return m_c1y;
  }

  CFreal getc1z() const
  {
    return m_c1z;
  }

  CFreal getc2x() const
  {
    return m_c2x;
  }

  CFreal getc2y() const
  {
    return m_c2y;
  }

  CFreal getc2z() const
  {
    return m_c2z;
  }

  CFreal getc3x() const
  {
    return m_c3x;
  }

  CFreal getc3y() const
  {
    return m_c3y;
  }

  CFreal getc3z() const
  {
    return m_c3z;
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
    return "LinearAdvSysTerm"; 
  }

protected:

  /// complete this list
  CFreal m_c0x;
  CFreal m_c0y;
  CFreal m_c0z;
  CFreal m_c1x;
  CFreal m_c1y;
  CFreal m_c1z;
  CFreal m_c2x;
  CFreal m_c2y;
  CFreal m_c2z;
  CFreal m_c3x;
  CFreal m_c3y;
  CFreal m_c3z;

}; // end of class LinearAdvSysTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdvSys_LinearAdvSysTerm_hh
