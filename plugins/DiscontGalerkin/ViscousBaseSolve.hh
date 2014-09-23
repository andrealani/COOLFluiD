#ifndef COOLFluiD_DiscontGalerkin_ViscousBaseSolve_hh
#define COOLFluiD_DiscontGalerkin_ViscousBaseSolve_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Framework/DataSocketSink.hh"
#include "DiscontGalerkin/DGElemTypeData.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"
#include "DiscontGalerkin/StdBaseSolve.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to assemble the system using Empty solver
 * @author Martin Holik
 */
class ViscousBaseSolve : public StdBaseSolve {

public: // functions

  /// Constructor
  explicit ViscousBaseSolve(const std::string& name);

  /// Destructor
  virtual ~ViscousBaseSolve();

  /// computes the K matrix
  void compute_Kmatrix2D(Framework::State m_state, std::vector < std::vector < RealMatrix > > *matrix);
  void compute_Kmatrix3D(Framework::State m_state, std::vector < std::vector < RealMatrix > > *matrix);
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

  CFreal m_Re;

}; // class ViscousBaseSolve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_ViscousBaseSolve_hh

