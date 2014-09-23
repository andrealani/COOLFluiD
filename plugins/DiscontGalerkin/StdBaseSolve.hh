#ifndef COOLFluiD_DiscontGalerkin_StdBaseSolve_hh
#define COOLFluiD_DiscontGalerkin_StdBaseSolve_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Framework/DataSocketSink.hh"
#include "DiscontGalerkin/DGElemTypeData.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to assemble the system using Empty solver
 * @author Martin Holik
 */
class StdBaseSolve : public DiscontGalerkinSolverCom {

public: // functions

  /// Constructor
  explicit StdBaseSolve(const std::string& name);

  /// Destructor
  virtual ~StdBaseSolve();

  /// computes matrix Pplus and Pminus for numerical flux
  CFuint compute_EigenValVec2D(Framework::State state, RealMatrix m_T, RealMatrix m_T1, RealMatrix *m_Pplus, RealMatrix *m_Pminus, RealMatrix *m_EigenVal, RealVector normal);
  CFuint compute_EigenValVec3D(Framework::State state, RealMatrix m_T, RealMatrix m_T1, RealMatrix *m_Pplus, RealMatrix *m_Pminus, RealMatrix *m_EigenVal, RealVector normal);

  /// computes matrix Pplus and Pminus for numerical flux
  //  void compute_EigenValVec_old_2D(Framework::State state, RealMatrix m_T, RealMatrix m_T1, RealMatrix *m_Pplus, RealMatrix *m_Pminus, RealMatrix *m_EigenVal, RealVector normal);

  /// computes the A matrix
  void compute_Amatrix2D(Framework::State m_state, std::vector < RealMatrix > *matrix);
  void compute_Amatrix3D(Framework::State m_state, std::vector < RealMatrix > *matrix);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

protected: // data

  ///temporary variable to store state;
  Framework::State *m_state;

}; // class StdBaseSolve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_StdBaseSolve_hh

