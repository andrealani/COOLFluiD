// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_BlendingP0Data_hh
#define COOLFluiD_FluxReconstructionMethod_BlendingP0Data_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/DataHandle.hh"
#include "Framework/State.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Reference element data of the cell-wise P0 order blending
 * (UseBlending = true), shared by the interior and boundary commands.
 *
 * Notation: c_j the cell-average weights of the solution points (sum 1),
 * divh_j,fp the divergence of the FR correction function of flux point fp at
 * solution point j, |J_j| the metric Jacobian at solution point j.
 *
 * - c_j: the cell-average weights of the element data (getCellAvgSolCoefs).
 * - W_fp = sum_j c_j divh_j,fp, over the solution points the FR correction of
 *   fp acts on. A face flux F_fp enters the cell average of the FR residual as
 *   -W_fp F_fp, so a term -sum_fp W_fp F_fp added at every solution point moves
 *   the cell average by the same amount.
 * - Cell average: the conserved average of the solution variables,
 *   sum_j c_j |J_j| U_j / sum_j c_j |J_j|, returned in the update variables
 *   through the physical data: computePhysicalData of the solution variable set,
 *   then computeStateFromPhysicalData of the update variable set.
 *
 * @author Rayan Dhib
 */
class BlendingP0Data {

public: // functions

  /// Constructor
  BlendingP0Data();

  /// Destructor
  ~BlendingP0Data();

  /**
   * Set up the weights and the variable sets.
   * @param frLocalData   element data of the cell type
   * @param corrFctDiv    divergence of the correction functions, [solution point][flux point]
   * @param flxSolDep     solution points each flux point depends on (the FR correction acts on these)
   * @param nbrEqs        number of equations
   * @param updateToSol   transformer from update to solution variables
   * @param solutionVar   solution variable set
   * @param updateVar     update variable set
   */
  void setup(FluxReconstructionElementData* frLocalData,
             const std::vector< std::vector< CFreal > >& corrFctDiv,
             Common::SafePtr< std::vector< std::vector< CFuint > > > flxSolDep,
             const CFuint nbrEqs,
             Common::SafePtr< Framework::VarSetTransformer > updateToSol,
             Common::SafePtr< Framework::ConvectiveVarSet > solutionVar,
             Common::SafePtr< Framework::ConvectiveVarSet > updateVar);

  /// Release the temporary states
  void unsetup();

  /// W_fp of cell flux point flxIdx
  CFreal getW(const CFuint flxIdx) const { return m_w[flxIdx]; }

  /// cell-average weight c_j of solution point solIdx
  CFreal getCellAvgCoef(const CFuint solIdx) const { return m_c[solIdx]; }

  /**
   * Conserved cell average of the states of one cell, in update variables.
   * @param states   states of the cell (update variables)
   * @param volumes  metric Jacobian |J_j| per state local ID
   * @param avg      output state
   */
  void computeCellAverage(const std::vector< Framework::State* >& states,
                          Framework::DataHandle< CFreal >& volumes,
                          Framework::State& avg);

private: // data

  /// number of solution points
  CFuint m_nbrSolPnts;

  /// number of equations
  CFuint m_nbrEqs;

  /// cell-average weights c_j
  std::vector< CFreal > m_c;

  /// W_fp per cell flux point
  std::vector< CFreal > m_w;

  /// transformer from update to solution variables (owned by the method data)
  Common::SafePtr< Framework::VarSetTransformer > m_updateToSol;

  /// solution variable set (owned by the method data)
  Common::SafePtr< Framework::ConvectiveVarSet > m_solutionVar;

  /// update variable set (owned by the method data)
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVar;

  /// physical data of the cell average
  RealVector m_pData;

  /// average of the solution variables
  Framework::State* m_solAvg;

}; // class BlendingP0Data

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BlendingP0Data_hh
