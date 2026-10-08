// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_SubcellBlendingQuadData_hh
#define COOLFluiD_FluxReconstructionMethod_SubcellBlendingQuadData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"
#include "Framework/Node.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"
#include "MathTools/RealVector.hh"
#include "Common/SafePtr.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Subcell grid of a quadrilateral element, shared by the four subcell P0
 * order blending commands (explicit and Jacobian, interior and boundary).
 *
 * One subcell per solution point. In each reference direction the subcell widths
 * are the quadrature weights of the 1D solution point distribution, so the
 * subcells tile the reference element and the P0 residual obeys the
 * same discrete conservation identity as the FR residual.
 *
 * Internal interfaces are stored as one flat list, the ones normal to ksi first
 * and the ones normal to eta after. Solution points are numbered
 * iSol = iKsi * nbrSolPnts1D + iEta, so the neighbours of a subcell in the two
 * reference directions are iSol +/- nbrSolPnts1D and iSol +/- 1.
 *
 * This is a plain data holder, not a command and not a configurable strategy.
 * Each command owns one instance.
 *
 * Subcell face vectors (computeCellNormals). The flux through an internal subcell face
 * is the Riemann flux in the direction of a face vector times its length. Along a line of
 * points in xi (fixed eta_j), let S_m be the xi-plane normal at solution point m (the
 * vector FR projects the flux on), S(-1), S(+1) the xi-plane normals at the two flux
 * points that end the line (the vectors of the FR face fluxes there), and
 * dS_k = sum_m D_km S_m + gL'_k (S(-1) - sum_m l_m(-1) S_m) + gR'_k (S(+1) - sum_m l_m(+1) S_m)
 * the derivative FR computes for S at point k: D_km is the derivative of the Lagrange
 * polynomial l_m of point m at point k, gL'_k and gR'_k the derivatives of the left and
 * right correction functions at point k. The face vector between the subcells of points
 * i and i+1 is n_(i+1/2) = S(-1) + sum_(k<=i) w_k dS_k (same along eta). The sum over the
 * whole line lands on S(+1), and the subcell residual of a uniform flux equals the FR
 * residual of that flux at the solution point, so the subcells keep a uniform flow exactly
 * whenever FR does. This is the telescoping of Hennemann et al. 2021 (appendix B.4, on
 * Gauss-Lobatto points, where the correction terms vanish), applied to the FR derivative
 * on Gauss-Legendre points. It needs the correction function (setCorrectionFunction);
 * without it the plane normal sampled at each face is used.
 *
 * Linear reconstruction (setReconstruction): along each reference line the subcell of
 * solution point i carries u(x) = u_i + s_i (x - x_i), x the reference coordinate of the
 * line and u the stored (update) variables. The slope s_i is a limited combination of the
 * secants to the two line neighbours. The end points of a line use as outer neighbour a
 * sample at the element face (x = -1 or +1): the trace of the neighbour cell at the face
 * flux point, or the boundary face value. The flux at a subcell face is the Riemann flux of
 * the two reconstructed states. The profile is exact for fields linear in the reference
 * coordinates, so for linear fields on affine cells, skewed ones included. The slopes and
 * the reconstructed values can be computed in other variables than the stored ones
 * (setReconstructionVars, e.g. primitive variables at high Mach). If a
 * reconstructed face state of a point is not admissible (update variable set isValid),
 * the slope of that point is set to zero.
 *
 * @author Rayan Dhib
 */
class SubcellBlendingQuadData {

public: // functions

  /// Constructor
  SubcellBlendingQuadData();

  /// Destructor
  ~SubcellBlendingQuadData();

  /// Build the subcell grid from the element data. Throws if the element is not a quad.
  void setup(FluxReconstructionElementData* frData, const CFuint dim, const CFuint nbrEqs);

  /// Quadrature weights of the 1D solution point distribution, used as subcell widths.
  /// Throws if a weight is not positive or if the weights do not sum to the reference length 2.
  static std::vector< CFreal > computeSubcellWidths1D(FluxReconstructionElementData* frData);

  /// number of solution points in one reference direction
  CFuint getNbrSolPnts1D() const { return m_nbrSolPnts1D; }

  /// subcell widths in one reference direction
  const std::vector< CFreal >& getWidths1D() const { return m_widths1D; }

  /// solution point of the subcell that owns the given cell flux point
  CFuint getClosestSol(const CFuint flxIdx) const { return (*m_closestSolToFlx)[flxIdx]; }

  /// width of that subcell in the direction normal to the face of the given flux point
  CFreal getFaceSubcellWidth(const CFuint flxIdx) const;

  /// Jacobian-scaled normals at the internal interfaces of the given cell, telescoped from
  /// the FR metric once setCorrectionFunction was called, sampled at the faces otherwise
  /// (see the class comment). Must be called once per cell before computeSubcellRes or
  /// addSubcellResDelta.
  void computeCellNormals(Framework::GeometricEntity* cell);

  /// Read the 1D correction function derivatives gL', gR' at the solution points from the
  /// divergence of the FR correction functions, corrFctDiv[solution point][cell flux point],
  /// and turn the telescoped normals on. Prints a warning and keeps the sampled normals if
  /// corrFctDiv does not have the tensor-product form expected on quads.
  void setCorrectionFunction(const std::vector< std::vector< CFreal > >& corrFctDiv);

  /// Turn the linear reconstruction on. limiter: "VanAlbada", "Minmod" or "None" (average
  /// of the two secants, for verification only). limiterEps: relative size of the van Albada
  /// smoothing (see limitSlope). updateVarSet gives the admissibility test.
  void setReconstruction(const std::string& limiter, const CFreal limiterEps,
                         Common::SafePtr< Framework::ConvectiveVarSet > updateVarSet);

  /// Reconstruct in other variables than the update ones: toRec maps update to
  /// reconstruction variables, fromRec back. Without this call the update variables are used.
  void setReconstructionVars(Common::SafePtr< Framework::VarSetTransformer > toRec,
                             Common::SafePtr< Framework::VarSetTransformer > fromRec);

  /// true if the linear reconstruction is on
  bool isLinear() const { return m_linear; }

  /// Limiter freezing. buffer: one limiter factor per (state local ID, direction, equation),
  /// valid: one flag per (state local ID, direction). A point whose flag is set uses the
  /// stored factor phi, slope = phi (a + b)/2 with a, b its two secants, instead of the live
  /// limiter. With record = true the live factor of every point not yet frozen is stored and
  /// flagged. Null buffers: live limiter everywhere.
  void setLimiterFreeze(CFreal* buffer, CFuint* valid, const bool record)
  {
    m_frozenPhi = buffer; m_frozenValid = valid; m_recordPhi = record;
  }

  /// Use the reconstruction in the current cell (true) or its first-order subcells (false).
  /// Only read when the reconstruction is on.
  void setCellLinear(const bool cellLinear) { m_cellLinear = cellLinear; }

  /// number of flux points of the element
  CFuint getNbrFlxPnts() const { return m_closestSolToFlx->size(); }

  /// Slopes of all solution points of the current cell along both reference directions.
  /// cellSamples: element face samples of the cell, layout [nbrEqs*flxIdx + iEq] over all
  /// cell flux points. The pointer is kept for addSubcellResDelta.
  /// Must be called before computeSubcellRes when the reconstruction is on.
  void computeCellSlopes(const std::vector< Framework::State* >& states, const CFreal* cellSamples);

  /// true if the last slope computed for this solution point along dir (KSI or ETA) kept the
  /// reconstruction, false if it was set to zero because a reconstructed state was not
  /// admissible; read right after computeCellSlopes
  bool isReconstructed(const CFuint dir, const CFuint sol) const { return m_reconstructed[dir][sol]; }

  /// Reconstructed state of the solution point closest to the given cell flux point, at
  /// the element face of that flux point. sample is the value at the face used as outer
  /// neighbour (the trace of the other cell). Used at element faces, where the cell loop
  /// slopes are not available yet; gives the same slope as computeCellSlopes.
  void reconstructAtElementFace(const std::vector< Framework::State* >& states,
                                const CFuint flxIdx, RealVector& sample,
                                Framework::State& result);

  /// Fill the per-solution-point subcell residual buffer from all internal interfaces,
  /// already scaled by alpha, and keep the interface fluxes for addSubcellResDelta.
  void computeSubcellRes(const CFreal alpha,
                         const std::vector< Framework::State* >& states,
                         RiemannFlux& riemannFlux);

  /// buffer filled by computeSubcellRes, indexed by solution point
  const std::vector< RealVector >& getSubcellRes() const { return m_subcellRes; }

  /// Apply the change in the subcell residual caused by perturbing one solution point,
  /// on top of the buffer of computeSubcellRes. Only the interfaces touching the
  /// perturbed point are recomputed, at most two per reference direction; with the
  /// reconstruction on, all interfaces of its two lines (frozen face samples). The residual
  /// uses the flat layout res[nbrEqs*iSol + iEq]. Both endpoints of those interfaces
  /// share a reference line with the perturbed point, so nothing outside the region
  /// filled by the parent perturbed volume term is written.
  void addSubcellResDelta(const CFreal alpha, const CFuint pertSol,
                          const std::vector< Framework::State* >& states,
                          RiemannFlux& riemannFlux,
                          RealVector& res);

private: // functions

  /// Riemann flux at an internal interface between the two adjacent solution point
  /// states, scaled by the size of the metric normal, written into the given vector.
  /// With the reconstruction on, the states are reconstructed with the given slopes
  /// (indexed by solution point) along the direction of the interface.
  void computeIntfFlux(const CFuint iIntf,
                       const std::vector< Framework::State* >& states,
                       const std::vector< RealVector >& slopes,
                       RiemannFlux& riemannFlux,
                       RealVector& result);

  /// recompute the flux of one internal interface with the given slopes and add the
  /// change with respect to the stored unperturbed flux to res (layout of addSubcellResDelta)
  void addIntfFluxDelta(const CFreal alpha, const CFuint iIntf,
                        const std::vector< Framework::State* >& states,
                        const std::vector< RealVector >& slopes,
                        RiemannFlux& riemannFlux,
                        RealVector& res);

  /// solution point number of point i (along dir) on the line with transverse index j
  CFuint lineSol(const CFuint dir, const CFuint j, const CFuint i) const
  {
    return (dir == KSI) ? i*m_nbrSolPnts1D + j : j*m_nbrSolPnts1D + i;
  }

  /// Limited slope of point i on the line (dir, j). lowSample/highSample: values at x = -1
  /// and x = +1, only read for the end points. The slope is zero if a reconstructed state
  /// at one of the two subcell faces of the point is not admissible; returns false then.
  bool computePointSlope(const CFuint dir, const CFuint j, const CFuint i,
                         const std::vector< Framework::State* >& states,
                         const CFreal* lowSample, const CFreal* highSample,
                         RealVector& slope);

  /// update variables to reconstruction variables (copy without transformers)
  void toRecVars(const RealVector& update, RealVector& rec);

  /// reconstruction variables to update variables (copy without transformers)
  void fromRecVars(const RealVector& rec, RealVector& update);

  /// Slope limiter applied to the two secants a and b. eps2 is the square of the van
  /// Albada smoothing size: with secants much smaller than it the limiter returns their
  /// average (a smooth extremum is not clipped), much larger ones are limited as usual
  /// (about zero for secants of opposite sign). The van Albada form is continuous and
  /// smooth near smooth extrema, which Newton needs to converge; Minmod is not.
  CFreal limitSlope(const CFreal a, const CFreal b, const CFreal eps2) const;

  /// Slopes of all points of the line (dir, j) into lineSlopes (indexed by solution
  /// point), with the samples of m_cellSamples.
  void computeLineSlopes(const CFuint dir, const CFuint j,
                         const std::vector< Framework::State* >& states,
                         std::vector< RealVector >& lineSlopes);

  /// Internal interface normals of the current cell telescoped from the plane normals at the
  /// solution points and at the element faces (see the class comment).
  void telescopeNormals(Framework::GeometricEntity* cell);

private: // data

  /// dimensionality of the physical model
  CFuint m_dim;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// number of solution points in one reference direction
  CFuint m_nbrSolPnts1D;

  /// subcell widths in one reference direction
  std::vector< CFreal > m_widths1D;

  /// solution point on the low side of each internal interface
  std::vector< CFuint > m_intfSolL;

  /// solution point on the high side of each internal interface
  std::vector< CFuint > m_intfSolR;

  /// subcell width of the low side of each internal interface
  std::vector< CFreal > m_intfWidthL;

  /// subcell width of the high side of each internal interface
  std::vector< CFreal > m_intfWidthR;

  /// mapped coordinates of each internal interface
  std::vector< RealVector > m_intfCoords;

  /// reference direction normal to each internal interface, KSI or ETA
  std::vector< CFuint > m_intfPlaneIdx;

  /// Jacobian-scaled normals at the internal interfaces of the current cell
  std::vector< RealVector > m_intfNormals;

  /// unperturbed flux at each internal interface of the current cell
  std::vector< RealVector > m_intfFlux;

  /// internal interfaces adjacent to each solution point
  std::vector< std::vector< CFuint > > m_intfOfSol;

  /// subcell residual per solution point, scaled by alpha
  std::vector< RealVector > m_subcellRes;

  /// perturbed fluxes of the interfaces adjacent to the perturbed solution point
  std::vector< RealVector > m_intfFluxPert;

  /// difference between the perturbed and unperturbed flux at one interface
  RealVector m_intfFluxDiff;

  /// solution point closest to each cell flux point
  Common::SafePtr< std::vector< CFuint > > m_closestSolToFlx;

  /// reference direction normal to the face of each cell flux point
  Common::SafePtr< std::vector< CFuint > > m_flxPntFlxDim;

  /// unit normal at the interface currently being evaluated
  RealVector m_unitNormal;

  /// true once setCorrectionFunction succeeded: the internal normals are telescoped
  bool m_telescope;

  /// 1D Lagrange derivative matrix at the solution points: [k][m] = derivative of the
  /// Lagrange polynomial of point m at point k
  std::vector< std::vector< CFreal > > m_derivMat1D;

  /// 1D Lagrange polynomials of the solution points at the line ends: [0][m] = l_m(-1),
  /// [1][m] = l_m(+1)
  std::vector< std::vector< CFreal > > m_lagrangeAtEnds;

  /// derivatives of the left and right correction functions at the 1D solution points
  std::vector< CFreal > m_corrDerivL;
  std::vector< CFreal > m_corrDerivR;

  /// mapped coordinates and plane index of the plane normals needed by the telescoping, in
  /// blocks: xi-plane at the solution points, eta-plane at the solution points, then per
  /// line xi-plane at xi = -1, eta-plane at eta = -1, xi-plane at xi = +1, eta-plane at eta = +1
  std::vector< RealVector > m_metricCoords;
  std::vector< CFuint > m_metricPlaneIdx;

  /// running sum of a telescoped face normal, the FR derivative of the plane normal at a
  /// point, and the two differences between face and extrapolated plane normals of a line
  RealVector m_faceSum;
  RealVector m_metricDeriv;
  RealVector m_endJumpL;
  RealVector m_endJumpR;

  /// true if the linear reconstruction is on
  bool m_linear;

  /// true if the current cell uses the reconstruction (see setCellLinear)
  bool m_cellLinear;

  /// frozen limiter factors and their flags (see setLimiterFreeze), null when not frozen
  CFreal* m_frozenPhi;
  CFuint* m_frozenValid;

  /// store the live limiter factors of points not frozen yet
  bool m_recordPhi;

  /// limiter: 0 Minmod, 1 VanAlbada, 2 None
  CFuint m_limiter;

  /// van Albada smoothing size relative to the largest of the three values of the stencil
  CFreal m_limiterEps;

  /// update variable set, for the admissibility test of reconstructed states
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;

  /// 1D solution point coordinates
  std::vector< CFreal > m_solPnts1D;

  /// 1D subcell boundaries: -1, -1 + w_0, ..., +1
  std::vector< CFreal > m_bnds1D;

  /// cell flux point at the end of each line: [dir][transverse index][0 low, 1 high]
  std::vector< std::vector< std::vector< CFuint > > > m_lineFlx;

  /// slopes of the current cell, [dir][solution point]
  std::vector< std::vector< RealVector > > m_slopes;

  /// false where the last slope of a line was set to zero by the admissibility test, [dir][solution point]
  std::vector< std::vector< bool > > m_reconstructed;

  /// slopes of a perturbed line, indexed by solution point
  std::vector< RealVector > m_pertSlopes;

  /// element face samples of the current cell (see computeCellSlopes)
  const CFreal* m_cellSamples;

  /// reconstructed left and right states at an internal interface
  Framework::State* m_recStateL;
  Framework::State* m_recStateR;

  /// coordinates given to the reconstructed states
  Framework::Node* m_recNodeL;
  Framework::Node* m_recNodeR;

  /// scratch values for the admissibility test
  RealVector m_testState;

  /// transformers between update and reconstruction variables (null: update variables)
  Common::SafePtr< Framework::VarSetTransformer > m_toRec;
  Common::SafePtr< Framework::VarSetTransformer > m_fromRec;

  /// scratch state given to the transformers
  Framework::State* m_transInState;

  /// reconstruction variables of the previous, current and next point of a stencil, and
  /// of a reconstructed face value
  RealVector m_recPrev;
  RealVector m_recCur;
  RealVector m_recNext;
  RealVector m_recFace;

  /// scratch slope of one point
  RealVector m_pointSlope;

}; // class SubcellBlendingQuadData

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_SubcellBlendingQuadData_hh
