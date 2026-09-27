// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_SubcellBlendingQuadData_hh
#define COOLFluiD_FluxReconstructionMethod_SubcellBlendingQuadData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
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

  /// Jacobian-scaled normals at the internal interfaces of the given cell, corrected so
  /// that every subcell is closed (see closeSubcells).
  /// Must be called once per cell before computeSubcellRes or addSubcellResDelta.
  void computeCellNormals(Framework::GeometricEntity* cell);

  /// Fill the per-solution-point subcell residual buffer from all internal interfaces,
  /// already scaled by alpha, and keep the interface fluxes for addSubcellResDelta.
  void computeSubcellRes(const CFreal alpha,
                         const std::vector< Framework::State* >& states,
                         RiemannFlux& riemannFlux);

  /// buffer filled by computeSubcellRes, indexed by solution point
  const std::vector< RealVector >& getSubcellRes() const { return m_subcellRes; }

  /// Apply the change in the subcell residual caused by perturbing one solution point,
  /// on top of the buffer of computeSubcellRes. Only the interfaces touching the
  /// perturbed point are recomputed, at most two per reference direction. The residual
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
  void computeIntfFlux(const CFuint iIntf,
                       const std::vector< Framework::State* >& states,
                       RiemannFlux& riemannFlux,
                       RealVector& result);

  /// Correct the internal interface normals of the current cell so that every subcell
  /// is closed: for each subcell the area vectors (normal x transverse width) of its four
  /// faces sum to zero, so a constant flux gives a zero subcell residual. The exterior
  /// faces (element flux points, shared with the neighbour element) are left untouched;
  /// the internal area vectors get the smallest change that closes all subcells.
  void closeSubcells(Framework::GeometricEntity* cell);

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

  /// false if the face flux points do not sit at the transverse solution point
  /// coordinates: then the subcell faces are not defined by the flux points and
  /// the internal normals are left as sampled
  bool m_closeSubcells;

  /// true once the warning about an element whose exterior faces do not close was printed
  bool m_warnedOpenElement;

  /// transverse subcell width of each internal interface
  std::vector< CFreal > m_intfTransWidth;

  /// mapped coordinates of the exterior subcell faces, 4 x nbrSolPnts1D:
  /// ksi = -1, ksi = +1, eta = -1, eta = +1, each ordered along the face
  std::vector< RealVector > m_extCoords;

  /// reference direction normal to each exterior subcell face, KSI or ETA
  std::vector< CFuint > m_extPlaneIdx;

  /// cosine eigenvectors of the 1D path graph Laplacian, [mode][point]
  std::vector< std::vector< CFreal > > m_cosine;

  /// eigenvalues of the 1D path graph Laplacian, per mode
  std::vector< CFreal > m_eigenvalue;

  /// closure defect of each subcell (sum of its outward area vectors)
  std::vector< RealVector > m_closureDefect;

  /// potential per subcell whose differences give the area vector corrections
  std::vector< RealVector > m_closurePotential;

  /// modal coefficient of the defect, scratch vector of size dim
  RealVector m_closureCoef;

}; // class SubcellBlendingQuadData

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_SubcellBlendingQuadData_hh
