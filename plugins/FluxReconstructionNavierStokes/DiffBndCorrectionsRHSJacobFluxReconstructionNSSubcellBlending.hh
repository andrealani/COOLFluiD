// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/DiffBndCorrectionsRHSJacobFluxReconstructionNS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Diffusive boundary command for NS that goes with
 * DiffRHSJacobFluxReconstructionNSSubcellBlending: for a boundary face of a cell
 * with blending coefficient alpha > 0, the face term of the cell becomes
 *   (1 - alpha) x FR correction + alpha x subcell boundary flux,
 * i.e. the parent's correction plus alpha sum_f F_f (c_f/w_f delta(closest(f)) - divh(., f)),
 * F_f the common boundary flux of the parent (ghost state, compact BR2 gradient)
 * times the face Jacobian size, c_f/w_f as in the interior command. The parent
 * calls computeCorrection for the residual and in its perturbation loop, so its
 * Jacobian includes the difference.
 *
 * Selected with UseSubcellDiffusion = true and SpaceRHSJacobCom = RHSJacobNS.
 *
 * @author Rayan Dhib
 */
class DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending : public DiffBndCorrectionsRHSJacobFluxReconstructionNS {

public: // functions

  /// Constructor
  explicit DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending() {}

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /// face term of the parent, (1 - alpha) x FR correction + alpha x subcell boundary flux
  /// in an alpha > 0 cell
  virtual void computeCorrection(std::vector< RealVector >& corrections);

protected: // data

  /// socket holding the blending coefficient of each state (one value per cell)
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// solution point closest to each cell flux point
  Common::SafePtr< std::vector< CFuint > > m_closestSolToFlx;

  /// c_f/w_f of each cell flux point
  std::vector< CFreal > m_faceFactors;

}; // class DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending_hh
