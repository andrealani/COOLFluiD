// Copyright (C) 2026 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_GammaAlphaTransitionCriterion_hh
#define COOLFluiD_FluxReconstructionMethod_GammaAlphaTransitionCriterion_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/DataHandle.hh"
#include "Framework/State.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Transition flag of the Gamma-Alpha model at the flux points of a boundary
 * face, set on the boundary condition for its wall ghost states. At every flux
 * point the flag is
 *
 *   tau / sqrt(rho mu) <= alpha
 *
 * with rho the density and mu the laminar dynamic viscosity of the state
 * extrapolated to the flux point, alpha the transition variable of that state
 * (component 5 + dim) and tau the wall shear built from the velocity gradient:
 * mu times the normal derivative of the tangential velocity, in 3D the norm
 * over the two tangent directions. The velocity gradient is the gradient of the
 * gradient variables of the interior cell corrected with this face only, at the
 * flux points: the volume term of the cell plus the lifting (g_b - g^D_f) grad h_f
 * of this face, with g^D_f the gradient variables of the cell extrapolated to
 * the flux points, g_b the boundary value the boundary condition lifts them to
 * (computeBndGradVars) and h_f the correction function of the face, divided by
 * the Jacobian determinant. g_b needs ghost states, which are computed first
 * with every flag false.
 *
 * Used by the convective boundary commands of the Gamma-Alpha model, which
 * differ only by their base class, once per face per residual evaluation; the
 * flags stay on the boundary condition for every later evaluation of the face.
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class GammaAlphaTransitionCriterion {

public: // functions

  /// Constructor
  GammaAlphaTransitionCriterion();

  /// Destructor
  ~GammaAlphaTransitionCriterion() {}

  /**
   * Sizes the scratch and stores the element data.
   * @param frData          element data of the cells
   * @param nbrEqs          number of equations
   * @param dim             dimension
   * @param updateVarSet    update variable set, for the pressure and temperature of a state
   * @param diffusiveVarSet diffusive variable set, for the gradient variables, the density and the viscosity
   */
  void setup(FluxReconstructionElementData& frData,
             const CFuint nbrEqs,
             const CFuint dim,
             Common::SafePtr< Framework::ConvectiveVarSet > updateVarSet,
             Common::SafePtr< Framework::DiffusiveVarSet > diffusiveVarSet);

  /**
   * Sets the flags of the current face on the boundary condition bc, whose face
   * must be set. The ghost states are overwritten.
   * @param cellStates              states of the interior cell
   * @param cellStatesFlxPnt        states extrapolated to the flux points of the face
   * @param ghostStates             ghost states at the flux points
   * @param unitNormals             unit normals at the flux points
   * @param flxPntCoords            coordinates of the flux points
   * @param faceJacobVecSizeFlxPnts face Jacobian vector sizes at the flux points
   * @param faceFlxPntConn          cell flux point index of every flux point of the face
   * @param nbrFaceFlxPnts          number of flux points of the face
   * @param corrFctDiv              divergence of the correction functions at the solution points
   * @param solPntNormals           metric of every state (dim x dim, row major)
   * @param volumes                 Jacobian determinant of every state
   */
  void setTransitionFlags(BCStateComputer& bc,
                          const std::vector< Framework::State* >& cellStates,
                          const std::vector< Framework::State* >& cellStatesFlxPnt,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          const std::vector< CFreal >& faceJacobVecSizeFlxPnts,
                          const std::vector< CFuint >& faceFlxPntConn,
                          const CFuint nbrFaceFlxPnts,
                          const std::vector< std::vector< CFreal > >& corrFctDiv,
                          const Framework::DataHandle< CFreal >& solPntNormals,
                          const Framework::DataHandle< CFreal >& volumes);

private: // functions

  /**
   * Wall shear at one flux point: mu times the normal derivative of the
   * tangential velocity, from the gradients of the velocity components.
   * @param normal   unit normal
   * @param grads    gradients of the gradient variables at the flux point
   * @param mu       laminar dynamic viscosity
   */
  CFreal wallShear(const RealVector& normal,
                   const std::vector< RealVector* >& grads,
                   const CFreal mu) const;

private: // data

  /// number of equations
  CFuint m_nbrEqs;

  /// dimension
  CFuint m_dim;

  /// number of solution points of a cell
  CFuint m_nbrSolPnts;

  /// number of solution points a solution point depends on
  CFuint m_nbrSolSolDep;

  /// solution points every solution point depends on
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solSolDep;

  /// derivatives of the solution polynomials at the solution points
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPnts;

  /// solution points every flux point depends on
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_flxSolDep;

  /// solution polynomials at the flux points
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;

  /// update variable set
  Common::SafePtr< Physics::NavierStokes::EulerVarSet > m_eulerVarSet;

  /// diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;

  /// the diffusive variable set as a Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

  /// physical data of a state
  RealVector m_pData;

  /// states of the cell, for setGradientVars
  std::vector< RealVector* > m_gradVarStatePtrs;

  /// gradient variables at the solution points
  RealMatrix m_gradVarsSolPnts;

  /// gradient variables extrapolated to the flux points, storage
  std::vector< RealVector > m_flxPntGradVarsStore;

  /// gradient variables extrapolated to the flux points, pointers into the storage
  std::vector< RealVector* > m_flxPntGradVars;

  /// boundary values of the gradient variables, storage
  std::vector< RealVector > m_bndGradVarsStore;

  /// boundary values of the gradient variables, pointers into the storage
  std::vector< RealVector* > m_bndGradVars;

  /// metric of the cell at the solution points, per direction
  std::vector< std::vector< RealVector > > m_cellFluxProjVects;

  /// scratch of the gradient kernels
  RealVector m_projectedCorr;

  /// gradient of the gradient variables corrected with the face, at the solution points
  std::vector< std::vector< RealVector > > m_gradsSolPnts;

  /// the same gradient at the flux points, storage
  std::vector< std::vector< RealVector > > m_gradsFlxPntStore;

  /// the same gradient at the flux points, pointers into the storage
  std::vector< std::vector< RealVector* > > m_gradsFlxPnt;

}; // class GammaAlphaTransitionCriterion

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_GammaAlphaTransitionCriterion_hh
