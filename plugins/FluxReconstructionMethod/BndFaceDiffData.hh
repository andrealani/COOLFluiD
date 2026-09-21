// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_BndFaceDiffData_hh
#define COOLFluiD_FluxReconstructionMethod_BndFaceDiffData_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "Common/COOLFluiD.hh"
#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Data of the diffusive boundary flux at the flux points of one boundary face,
 * as computed by DiffBndCorrectionsRHSFluxReconstruction::computeBndFaceDiffData.
 * The first index of every member is the flux point.
 *
 * @author Rayan Dhib
 */
struct BndFaceDiffData
{
  /// interior states extrapolated to the flux points
  std::vector< RealVector > intStates;

  /// ghost states
  std::vector< RealVector > ghostStates;

  /// boundary states the diffusive flux is evaluated at
  std::vector< RealVector > bndStates;

  /// gradient variables extrapolated to the flux points
  std::vector< RealVector > gradVarsFlxPnt;

  /// boundary values of the gradient variables
  std::vector< RealVector > bndGradVars;

  /// boundary gradients [iFlx][iEq]
  std::vector< std::vector< RealVector > > bndGrads;

  /// diffusive fluxes, without convective and artificial viscosity terms
  std::vector< RealVector > diffFluxes;

  /// unit normals, pointing out of the fluid
  std::vector< RealVector > unitNormals;

  /// coordinates of the flux points
  std::vector< RealVector > coords;

  /// face Jacobian vector sizes
  std::vector< CFreal > faceJacobVecAbsSizes;

  /// face integration coefficients
  std::vector< CFreal > faceIntegrationCoefs;
};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BndFaceDiffData_hh
