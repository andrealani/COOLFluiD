#ifndef COOLFluiD_DiscontGalerkin_ViscousSolveCells_hh
#define COOLFluiD_DiscontGalerkin_ViscousSolveCells_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Framework/DataSocketSink.hh"

#include "DiscontGalerkin/DGElemTypeData.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"
#include "DiscontGalerkin/ViscousBaseSolve.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to assemble the system using Empty solver
 * @author Martin Holik
 */
class ViscousSolveCells : public ViscousBaseSolve{
public: // functions

  /// Constructor
  explicit ViscousSolveCells(const std::string& name);

  /// Destructor
  virtual ~ViscousSolveCells();

  /// Execute processing actions
  void execute();

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

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets();

  /**
   * Function for setting time step
   * @return a value of time step
   */
  CFreal setTimeStep(CFreal tau);

  /// map of LSSMatrix accumulators, one for each cell type
  Common::CFMap<CFuint,DGElemTypeData> m_mapElemData;
  /// socket for Rhs
  Framework::DataSocketSink<CFreal> socket_rhs;
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State*, Framework::GLOBAL >  socket_states;
  /// the socket to the data handle of the old state's
  Framework::DataSocketSource < Framework::State* >  socket_old_states;

private: // data

  /// handle for the InnerCells trs
  Common::SafePtr<Framework::TopologicalRegionSet> m_cells;

  /*                   /  invJacobi[0]  invJacobi[1]  \
  //     invJacobi =  |                                |
  //                   \  invJacobi[2]  invJacobi[3]  /
  */
  /// inverse Jacobi matrix
  CFreal invJacobi[4];

  /// determinant of Jacobi matrix
  CFreal detJacobi;

  /// matrix of linearization of viscous terms K[s][k](x,x)
  std::vector < std::vector < RealMatrix > > m_kMatrix;

  /// jacobi matrix of inviscid terms A[s](x,x)
  std::vector < RealMatrix > m_aMatrix;

  ///temporary variable to store old state;
  Framework::State *m_oldState;

  CFreal m_Alpha;
  CFreal m_MaxCFL;

}; // class ViscousSolveCells

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_ViscousSolveCells_hh

