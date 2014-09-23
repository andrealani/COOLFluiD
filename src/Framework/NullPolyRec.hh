// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullPolyRec_hh
#define COOLFluiD_Framework_NullPolyRec_hh

//////////////////////////////////////////////////////////////////////////////

#include "PolyReconstructor.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class implements a constant polynomial reconstructor for FVM
/// @author Andrea Lani
template < typename METHODDATA >
class NullPolyRec : public Framework::PolyReconstructor<METHODDATA> {
public:

  /// Constructor
  NullPolyRec(const std::string& name);

  /// Default destructor
  ~NullPolyRec();

  /// Compute the gradients
  void computeGradients();
  
protected:
  
  /// Allocate reconstruction data needed for the flux evaluation
  void allocateReconstructionData();
  
  /// Compute the flux in the current face
  void extrapolateImpl(Framework::GeometricEntity* const geo);

  /// Compute the flux in the current face
  void extrapolateImpl(Framework::GeometricEntity* const geo,
  	       CFuint iVar, CFuint leftOrRight);

  /// Compute the flux in the current face
  void baseExtrapolateImpl(Framework::GeometricEntity* const geo);

  /// Compute the flux in the current face
  void baseExtrapolateImpl(Framework::GeometricEntity* const geo,
  		   CFuint iVar, CFuint leftOrRight);

}; // end of class NullPolyRec

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

#include "NullPolyRec.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullPolyRec_hh
