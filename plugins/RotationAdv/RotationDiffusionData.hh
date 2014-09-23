// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdv_RotationDiffusionData_hh
#define COOLFluiD_Physics_RotationAdv_RotationDiffuionData_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics { 

    namespace RotationAdv { 
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the data for the NavierStokes physical model
 * 
 * @author Andrea Lani
 */
class RotationDiffusionData {

public:
  
  /**
   * Overloading of the assignment operator
   */
  const RotationDiffusionData& operator=(const RotationDiffusionData& data) 
  {
    nu = data.nu; // dynamic viscosity

    return *this;
  }
  
  /**
   * Set up private data
   */
  void setup()
  {
  }
  
  CFreal nu; // dynamic viscosity
  
};

//////////////////////////////////////////////////////////////////////////////
 
    } // namespace NavierStokes
 
  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokesData_hh
