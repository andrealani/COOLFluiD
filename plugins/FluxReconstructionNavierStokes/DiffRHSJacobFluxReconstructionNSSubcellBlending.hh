// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionNSSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionNSSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/DiffRHSJacobFluxReconstructionNS.hh"
#include "FluxReconstructionMethod/SubcellBlendingQuadData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

    class FluxReconstructionElementData;

//////////////////////////////////////////////////////////////////////////////

/**
 * Implicit diffusive RHS command for NS with subcell blending of the diffusive
 * term, for quads: in a cell with blending coefficient alpha > 0 (socket alpha,
 * written by the order blending)
 *   R_diff = (1 - alpha) R_BR2 + alpha R_sub,
 * the diffusive counterpart of the convective subcell blending, on the same
 * subcells (one per solution point, widths = Gauss-Legendre weights w) and the
 * same telescoped subcell normals (SubcellBlendingQuadData). Selected with
 * UseSubcellDiffusion = true next to UseSubcellBlending = true: the solver then
 * appends SubcellBlending to DiffSolveCom = DiffRHSJacobNS and to the diffusive
 * boundary command (DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending).
 *
 * R_sub is a finite-volume diffusion with two-point fluxes. At solution point i
 * of a reference line, R_sub,i = (F_(i+1/2) - F_(i-1/2)) / w_i, with F the
 * diffusive flux projected on the Jacobian-scaled subcell normal. At the element
 * faces F is the common face flux, which is itself blended,
 *   F* = (1 - alpha_F) F*_BR2 + alpha_F F_2p, alpha_F = max(alpha_L, alpha_R),
 * so that it stays single valued (as FaceFluxBlending for the convective flux).
 *
 * Two-point flux F_2p between two points L and R (neighbouring solution points
 * of a cell, or the solution points closest to an element face flux point in the
 * two cells): the diffusive flux of the var set at the face state (average of the
 * two update states) with, for each gradient variable u, the gradient
 *   g_f = gbar + (g2p - gbar.e) e,
 * e the unit vector from L to R, d their distance, gbar the average of the two
 * point gradients and g2p the derivative along e from the two values:
 * (u_R - u_L)/d for a linear variable, (exp(u_R - u_f) - exp(u_L - u_f))/d for
 * a logarithm (ln y_i, ln T, ln Tv of LogRhoivLogTTv), u_f the value at the face
 * state. The var set forms grad(y) = y_f grad(ln y) with the face value y_f, so the
 * logarithms get exactly the linear two-point difference (y_R - y_L)/d.
 *
 * Point gradient of a point of an alpha > 0 cell: the derivatives of u along the
 * two reference lines through the point, central differences of the neighbouring
 * points of the line (at an inner element face the face sample, the average of the
 * two solution points closest to the face flux point, placed at -1 or +1; one-sided
 * at a boundary face), mapped as in FR, grad u = (n_ksi du/dksi + n_eta du/deta)/|J|,
 * n the mapped coordinate plane normals and |J| the Jacobian determinant. Points
 * of alpha = 0 cells use the BR2 gradient.
 *
 * Implementation: computeDivDiscontFlx and computeDivDiscontFlxNeighb blend the
 * volume term with the internal subcell fluxes, computeInterfaceFlxCorrection
 * blends the common face flux, and computeCorrection gives the face term of an
 * alpha > 0 cell,
 *   (1 - alpha) sum_f F_f divh(., f) + alpha sum_f F_f c_f/w_f delta(closest(f)),
 * c_f = sum along the line of f of w_i divh(i, f): the subcell boundary flux of the
 * solution point closest to f instead of the FR correction, with the same integral.
 * The parent calls these functions for the residual and in its perturbation loops,
 * so the assembled Jacobian is that of the blended residual, except that gbar is
 * frozen in the perturbations.
 * The boundary faces are handled by DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending.
 *
 * @author Rayan Dhib
 */
class DiffRHSJacobFluxReconstructionNSSubcellBlending : public DiffRHSJacobFluxReconstructionNS {

public: // functions

  /// Constructor
  explicit DiffRHSJacobFluxReconstructionNSSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~DiffRHSJacobFluxReconstructionNSSubcellBlending() {}

  /// subcell data of the alpha > 0 cells, then residual and Jacobian of the parent
  virtual void execute();

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unset up private data
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

  /// For each cell flux point f, c_f/w_f: c_f = sum over the solution points i of
  /// the reference line of f of w_i divh(i, f), w_f the subcell width of the solution
  /// point closest to f along the face normal. Throws if c_f is not +-1 or if divh
  /// does not vanish outside the line (tensor-product correction functions).
  static void computeFaceSubcellFactors(FluxReconstructionElementData* frData,
                                        const std::vector< std::vector< CFreal > >& corrFctDiv,
                                        std::vector< CFreal >& factors);

protected: // functions

  /// common diffusive flux of the parent, blended with the two-point flux by alpha_F
  virtual void computeInterfaceFlxCorrection();

  /// face term of the parent, (1 - alpha) x FR correction + alpha x subcell boundary flux
  /// in an alpha > 0 cell (see the class comment)
  virtual void computeCorrection(CFuint side, std::vector< RealVector >& corrections);

  /// (1 - alpha) x parent volume term + alpha x internal subcell term of the cell of m_cellStates
  virtual void computeDivDiscontFlx(std::vector< RealVector >& residuals);

  /// (1 - alpha) x parent perturbed volume term + alpha x internal subcell term of the cell of side,
  /// updated for the perturbed point m_pertSol when side == m_pertSide
  virtual void computeDivDiscontFlxNeighb(RealVector& residuals, const CFuint side);

  /// face samples of the cells of every inner face with alpha > 0 on one side
  void computeFaceSamples();

  /// difference quotient of gradient variable iEq along reference direction dir at
  /// position i of line j of the cell in m_blendGradVars
  CFreal lineDerivative(const CFuint dir, const CFuint j, const CFuint i, const CFuint iEq, const CFuint cellID) const;

  /// for each cell with alpha > 0: subcell gradients, subcell normals, internal fluxes and R_sub
  void prepareAlphaCells();

  /// internal subcell fluxes and R_sub (without alpha) of a cell from its current states
  void computeCellSubcellRes(const CFuint cellID, const std::vector< Framework::State* >& states);

  /// two-point flux of internal interface iIntf of a cell, projected on its subcell normal
  void computeIntfFlux(const CFuint cellID, const std::vector< Framework::State* >& states,
                       const CFuint iIntf, RealVector& flux);

  /// two-point diffusive flux between states sL and sR, projected on normal (not unit),
  /// gradL and gradR the point gradients giving gbar
  void computeTwoPointFlux(Framework::State& sL, Framework::State& sR,
                           const std::vector< RealVector >& gradL, const std::vector< RealVector >& gradR,
                           const RealVector& normal, const CFreal radius, RealVector& flux);

  /// point gradient used for gbar: subcell gradient for a point of an alpha > 0 cell, BR2 otherwise
  const std::vector< RealVector >& getPointGrad(Framework::State& state, const std::vector< RealVector >& br2Grad) const;

protected: // data

  /// subcell data of one cell, cached per evaluation
  struct CellSubcellData {
    CFuint stamp;                          ///< evaluation of intfFlux and subRes
    std::vector< RealVector > normals;     ///< telescoped internal subcell normals (geometry, kept)
    std::vector< RealVector > intfFlux;    ///< unperturbed internal two-point fluxes
    std::vector< RealVector > subRes;      ///< unperturbed R_sub of the internal faces, without alpha
  };

  /// socket holding the blending coefficient of each state (one value per cell)
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// subcell grid (widths, telescoped normals)
  SubcellBlendingQuadData m_scData;

  /// number of solution points in one reference direction
  CFuint m_nbrSolPnts1D;

  /// 1D solution point coordinates and subcell widths
  std::vector< CFreal > m_solPnts1D;
  std::vector< CFreal > m_widths1D;

  /// solution point closest to each cell flux point
  Common::SafePtr< std::vector< CFuint > > m_closestSolToFlx;

  /// number of flux points of a cell
  CFuint m_nbrCellFlxPnts;

  /// cell flux point at the ends of each reference line, [dir][end][line], end 0 at -1, 1 at +1
  std::vector< std::vector< std::vector< CFuint > > > m_lineEndFlx;

  /// c_f/w_f of each cell flux point
  std::vector< CFreal > m_faceFactors;

  /// low and high solution point and their subcell widths of each internal interface
  std::vector< CFuint > m_intfSolL;
  std::vector< CFuint > m_intfSolR;
  std::vector< CFreal > m_intfWidthL;
  std::vector< CFreal > m_intfWidthR;

  /// internal interfaces of each solution point
  std::vector< std::vector< CFuint > > m_intfOfSol;

  /// true for the gradient variables that are logarithms
  std::vector< bool > m_isLogGradVar;

  /// cell of each local state
  std::vector< CFuint > m_cellOfState;

  /// face samples, [(cellID*m_nbrCellFlxPnts + flxIdx)*nbrEqs + iEq], and their evaluation
  std::vector< CFreal > m_faceSamples;
  std::vector< CFuint > m_sampleStamp;

  /// subcell gradients per local state [stateID][iEq], and their evaluation
  std::vector< std::vector< RealVector > > m_subGrads;
  std::vector< CFuint > m_subGradStamp;

  /// subcell data per cell
  std::vector< CellSubcellData > m_cellData;

  /// counter of the evaluations with alpha > 0 somewhere
  CFuint m_stamp;

  /// gradient variables of one cell, nbrEqs x nbrSolPnts, and of the two sides of a face
  RealMatrix m_blendGradVars;
  std::vector< RealMatrix > m_sideGradVars;

  /// mapped coordinate plane normals at the solution points of one cell, [dir][sol]
  std::vector< std::vector< RealVector > > m_blendProjVects;

  /// two-point flux work data: states, their gradient variables, gradients, unit vector, flux
  std::vector< RealVector* > m_twoStatePtrs;
  RealMatrix m_twoGradVars;
  std::vector< RealVector* > m_twoGrad;
  RealVector m_twoDir;
  RealVector m_twoFlux;
  RealVector m_twoFluxPert;

}; // class DiffRHSJacobFluxReconstructionNSSubcellBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionNSSubcellBlending_hh
