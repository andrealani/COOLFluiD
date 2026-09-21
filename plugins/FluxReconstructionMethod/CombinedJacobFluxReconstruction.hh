// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_CombinedJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_CombinedJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BlockAccumulator.hh"

#include "FluxReconstructionMethod/DiffRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Base class of the commands that compute a combined convective, diffusive and
 * artificial viscosity residual together with its Jacobian.
 *
 * Notation, the same for the fluxes and for the gradients. For one cell, F^D
 * is the discontinuous flux, the flux of the cell polynomial, F^D_f its value
 * extrapolated to the flux points of face f, F^I_f the interface flux of face f
 * (the Riemann flux, or the boundary flux) and h_f the correction function of
 * face f. The divergence of the corrected flux is
 *
 *   div F = div F^D + sum_f (F^I_f - F^D_f) div h_f
 *
 * The gradient variables g = g(U) of the diffusive variable set are treated the
 * same way: g^D is the polynomial through their values at the solution points,
 * g^D_f its value extrapolated to the flux points of face f, and g^I_f the
 * interface value of face f: the average of the two sides at an interior face,
 * the boundary value of the boundary condition at a boundary face. The
 * reconstructed gradient of the cell, corrected with all its faces, is
 *
 *   q = grad g^D + sum_f (g^I_f - g^D_f) grad h_f
 *
 * with grad h_f the unit normal of the face, scaled by the face Jacobian, times
 * div h_f. The compact gradient of one face f keeps the correction of that face
 * only, scaled by eta, and is taken at the flux points of that face,
 *
 *   q_f = grad g^D + eta (g^I_f - g^D_f) grad h_f
 *
 * so it depends on the two cells of the face alone. The interface diffusive
 * flux of face f is the diffusive flux at the average of the two extrapolated
 * states and at the average of the two sides' q_f; the volume flux of the cell
 * uses q. eta is the BR2 lifting
 * multiplier (option BR2Eta, 5 by default): it sets how strongly the jump of
 * the face enters the gradient the interface flux sees, and the usual BR2
 * condition wants it larger than the number of faces of the element. Both
 * gradients are built in mapped coordinates and divided by the Jacobian
 * determinant.
 *
 * The Jacobian is assembled with the chain rule. The partial derivatives of
 * the flux at the solution points and of the interface flux with respect to
 * the state and to the gradient are computed numerically. When one state is
 * perturbed, the changes of q and of q_f follow from the two expressions above,
 * which are linear in g^D and g^I, and the flux changes built from them are
 * sent through div F^D and the face corrections (F^I_f - F^D_f) div h_f.
 *
 * @author Rayan Dhib
 */
class CombinedJacobFluxReconstruction : public DiffRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit CombinedJacobFluxReconstruction(const std::string& name);

  /// Destructor
  virtual ~CombinedJacobFluxReconstruction() {}

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

protected: // functions

  /// @return true when the command has a physical diffusive flux to differentiate
  virtual bool hasPhysicalDiffusionJacobian() const
  {
    return true;
  }

  /// @return true when the command has an artificial viscosity flux to differentiate
  virtual bool hasArtificialViscosityJacobian() const
  {
    return false;
  }

  /// @return true when the artificial viscosity flux is also applied at the boundary faces
  virtual bool hasAVBoundaryFlux() const
  {
    return false;
  }

  /**
   * Set the artificial viscosity data of a cell whose faces are all boundary
   * faces, before its residual and Jacobian are computed.
   * @pre m_cells[LEFT] and m_states[LEFT] hold the cell
   */
  virtual void prepareIsolatedCellAV(const CFuint cellID)
  {
  }

  /**
   * compute the contribution of the diffusive face term to both Jacobians
   */
  virtual void computeBothJacobsDiffFaceTerm();

  /**
   * compute the contribution of the diffusive face term to one Jacobians
   */
  virtual void computeOneJacobDiffFaceTerm(const CFuint side);

  /**
   * Compute the derivative of the common face flux with respect to the states
   * extrapolated to the flux points, times the residual factor.
   */
  virtual void computeRiemannFluxJacobianNum(const CFreal resFactor) = 0;

  /**
   * Compute the derivative of the common face flux with respect to the average
   * face gradient, times the residual factor.
   */
  virtual void computeRiemannFluxToGradJacobianNum(const CFreal resFactor) = 0;

  /**
   * Compute, at every solution point of one cell, the derivative of the flux
   * with respect to the state (m_fluxJacobian) and to the gradients
   * (m_gradientFluxJacobian), both times minus the residual factor. The flux is
   * the diffusive flux minus the convective flux, the gradients are the
   * gradients corrected with all the faces of the cell.
   */
  void computeCellFluxJacobians(const CFuint side);

  /**
   * Assemble the Jacobian blocks of the current interior face.
   * @param ownedSide side whose rows are assembled, or -1 for both sides
   */
  void assembleFaceJacobian(const CFint ownedSide);

  /**
   * Put in m_cellGrads the change dq of the gradient corrected with all the
   * faces, for both cells, when variable m_pertVar at solution point m_pertSol
   * of cell m_pertSide is perturbed: the volume term grad dg^D of the perturbed
   * cell and the corrections (dg^I_f - dg^D_f) grad h_f of every face of that
   * cell, with dg^D the change of its gradient variables.
   * @return the inverse of the perturbation
   */
  CFreal computePertCellGradients();

  /**
   * Put in m_cellGradFlxPnt the change dq_f of the compact gradient of the
   * current face for the same perturbation: the volume term grad dg^D of the
   * perturbed cell and eta times the correction (dg^I_f - dg^D_f) grad h_f of
   * this face, taken at its flux points.
   */
  void computePertCompactFaceGradients(const CFuint nbrFaceFlxPnts, const CFreal invEps);

  /**
   * Add to the accumulator the derivative of the volume residual of one cell.
   * At every solution point the flux change is dF = F_U dU + F_q dq, with F_U
   * and F_q the partial derivatives of the flux with respect to the state and
   * to the gradient; it is sent through div dF^D and the face corrections
   * -dF^D_f div h_f of the cell. The interface flux part of the corrections is
   * added by addFaceFluxJacobian when the face is the current face.
   */
  void addCellVolumeJacobian(Framework::BlockAccumulator& acc, const CFuint destSide, const CFreal invEps);

  /**
   * Add to the face block accumulator the derivative of the face correction of
   * one cell. At every flux point of the current face the change of the
   * interface flux is dF^I_f = F^I_U dU^D_f + F^I_q dq_avg, with F^I_U and
   * F^I_q the partial derivatives of the interface flux, dU^D_f the state change
   * extrapolated to the flux point and dq_avg the average of the two compact
   * gradient changes; it is multiplied by the face Jacobian size and by div h_f
   * at every solution point of the cell.
   */
  void addFaceFluxJacobian(const CFuint destSide, const CFuint nbrFaceFlxPnts);

  /**
   * Compute the residual and assemble the Jacobian block of every owned cell
   * whose faces are all boundary faces: the face loop never visits them.
   */
  void computeCellsWithoutInnerFace();

  /**
   * Residual and Jacobian block of a cell whose faces are all boundary faces:
   * its volume residual with the gradient corrected with all its faces.
   * @pre m_cells[LEFT] and m_states[LEFT] hold the cell
   */
  void assembleIsolatedCellJacobian(const CFuint cellID);

  /**
   * Compute the derivatives of the artificial viscosity gradients for the
   * perturbation of variable iVar at solution point iSol of cell side: the
   * gradient corrected with all the faces (m_avGradDerivs), the compact gradient
   * of the current face (m_avFaceGradDerivs) and the derivative of the boundary
   * artificial viscosity residual (m_avBndResDerivs). The artificial viscosity is
   * frozen and the gradient variables are the conservative variables.
   * @param isolated true for a cell whose faces are all boundary faces
   */
  void computeAVGradientDerivatives(const CFuint side, const CFuint iSol, const CFuint iVar, const bool isolated);

  /**
   * Add to the flux change at solution point iSol of cell side, in direction
   * iDim, the artificial viscosity term: minus the residual factor, times the
   * viscosity eps of that solution point, times the change of the artificial
   * viscosity gradient projected on the mapped coordinate plane normal of the
   * direction.
   */
  void addAVCellFluxDerivative(RealVector& derivFlux, const CFuint side, const CFuint iSol, const CFuint iDim);

  /**
   * Add to the change of the common face flux at flux point iFlx the artificial
   * viscosity term: minus the residual factor, times the average viscosity of
   * the two sides, times the average change of the two compact artificial
   * viscosity gradients projected on the unit normal.
   */
  void addAVFaceFluxDerivative(RealVector& derivFlux, const CFuint iFlx);

  /**
   * Add to the derivative of the residual at solution point iSol of cell side
   * the derivative of the boundary artificial viscosity residual.
   */
  void addAVBndResidualDerivative(RealVector& derivRes, const CFuint side, const CFuint iSol);

protected: // data

  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;

  /// artificial viscosity at the solution points of both cells, [side][iSol]
  std::vector< std::vector< CFreal > > m_solEpsilons;

  /// artificial viscosity extrapolated to the flux points of the current face, [side][iFlx]
  std::vector< std::vector< CFreal > > m_epsilonLR;

  /// derivative of the flux with respect to the state, [side][iSol][iVar][iDim]
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_fluxJacobian;

  /// derivative of the flux with respect to the gradients, [side][iSol][iEq][iGradDim][iDim]
  std::vector< std::vector< std::vector< std::vector< std::vector< RealVector > > > > > m_gradientFluxJacobian;

  /// derivative of the common face flux with respect to the extrapolated states, [side][iFlx][iVar]
  std::vector< std::vector< std::vector< RealVector > > > m_riemannFluxJacobian;

  /// derivative of the common face flux with respect to the average face gradient, [iFlx][iEq][iDim]
  std::vector< std::vector< std::vector< RealVector > > > m_riemannFluxGradJacobian;

  /// common diffusive flux at the flux points of the current face
  std::vector< RealVector > m_flxPntRiemannFluxDiff;

  /// physical data of a state
  RealVector m_pData;

  /// flux at the current solution point in every direction, unperturbed
  std::vector< RealVector > m_unpertContFlx;

  /// flux at the current solution point in every direction, perturbed
  std::vector< RealVector > m_pertContFlx;

  /// derivative of the flux at one solution point in one direction
  RealVector m_derivContFlx;

  /// derivative of the flux at every solution point of a cell, [iSol][iDim]
  std::vector< std::vector< RealVector > > m_derivContFlxSolPnts;

  /// derivative of the residual at one solution point
  RealVector m_derivSolPntRes;

  /// derivative of a flux at one flux point
  RealVector m_derivFlxPntFlux;

  /// derivative of the gradient variables at the solution points of both cells, [side]
  std::vector< RealMatrix > m_derivGradVarsSolPnts;

  /// that derivative extrapolated to the flux points of the current face, [side](iEq,iFlx)
  std::vector< RealMatrix > m_derivGradVarsFlxPnt;

  /// local indexes of all the faces of a cell
  std::vector< CFuint > m_allFaceLocalIdxs;

  /// backup of the pointers to the gradients of both cells
  std::vector< std::vector< std::vector< RealVector >* > > m_cellGradsPtrsBackUp;

  /// backup of the physical gradient variables before the perturbation
  RealMatrix m_physGradVarsSolPntsBefore;

  /// derivative of the artificial viscosity gradient corrected with all the faces, [side][iSol][iEq]
  std::vector< std::vector< std::vector< RealVector > > > m_avGradDerivs;

  /// derivative of the artificial viscosity compact gradient of the current face, [side][iFlx][iEq]
  std::vector< std::vector< std::vector< RealVector > > > m_avFaceGradDerivs;

  /// pointers to m_avFaceGradDerivs, [side][iFlx][iEq]
  std::vector< std::vector< std::vector< RealVector* > > > m_avFaceGradDerivPtrs;

  /// derivative of the boundary artificial viscosity residual, [side] with size nbrSolPnts*nbrEqs
  std::vector< RealVector > m_avBndResDerivs;

  /// boundary artificial viscosity residual of a cell, unperturbed
  RealVector m_avBndRes;

  /// boundary artificial viscosity residual of a cell, perturbed
  RealVector m_pertAVBndRes;

}; // class CombinedJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_CombinedJacobFluxReconstruction_hh
