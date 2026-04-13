// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_BaseOrderBlendingMHD_hh
#define COOLFluiD_FluxReconstructionMethod_BaseOrderBlendingMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/BaseOrderBlending.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * MHD-aware order blending command. Extends BaseOrderBlending with the
 * "B2" monitored expression (magnetic pressure), accessed via the physical
 * data slots defined by MHDTerm — no raw state-index dependencies.
 *
 * All other expressions (rho, p, rho*p, p/rho, rho/p, velocity_magnitude)
 * are delegated to the base class.
 *
 * @author Rayan Dhib
 */
class BaseOrderBlendingMHD : public BaseOrderBlending {
public:

  explicit BaseOrderBlendingMHD(const std::string& name);

  virtual ~BaseOrderBlendingMHD();

protected:

  /// Handle "B2" via MHDTerm::BX/BY/BZ physical data slots; delegate otherwise.
  virtual void extractMonitoredField();

}; // class BaseOrderBlendingMHD

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BaseOrderBlendingMHD_hh
