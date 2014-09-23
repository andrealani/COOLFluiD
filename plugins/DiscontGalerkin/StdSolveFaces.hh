#ifndef COOLFluiD_DiscontGalerkin_StdSolveFaces_hh
#define COOLFluiD_DiscontGalerkin_StdSolveFaces_hh

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
 * This is a standard command to assemble the system using
 * discontinuous Galerkin solver by adding the contributions integrals
 * in the faces.
 *
 * @author Martin Holik
 */
class StdSolveFaces : public StdBaseSolve
{

public: // functions

  /// Constructor
  explicit StdSolveFaces(const std::string& name);

  /// Destructor
  virtual ~StdSolveFaces();

  /// Execute processing actions
  virtual void execute();

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

private :

  /// handle for the InnerCells trs
  Common::SafePtr<Framework::TopologicalRegionSet> m_cells;

  /// map of LSSMatrix accumulators, one for each cell type
  Common::CFMap<CFuint,DGElemTypeData> m_mapElemData;

  /// socket for Rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for integration indexes
  Framework::DataSocketSink< std::vector< CFuint > >
    socket_integrationIndex;

  /// socket for normals on face
  Framework::DataSocketSink< RealVector >
    socket_normals;

    /*                   /  invJacobi[0]  invJacobi[1]  \
    //     invJacobi =  |                                |
    //                   \  invJacobi[2]  invJacobi[3]  /
    */
  /// inverse Jacobi matrix
  //CFreal invJacobiLeft[4];
  //CFreal invJacobiRight[4];

  /// determinant of Jacobi matrix
  //CFreal detJacobiLeft;
  //CFreal detJacobiRight;

}; // class StdSolveFaces

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_StdSolveFaces_hh

