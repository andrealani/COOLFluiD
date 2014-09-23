#ifndef COOLFluiD_DiscontGalerkin_WallBC_hh
#define COOLFluiD_DiscontGalerkin_WallBC_hh

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
 * This is a standard command to assemble the system using
 * discontinuous Galerkin solver by adding the contributions integrals
 * in the faces.
 *
 * @author Martin Holik
 */
class WallBC : public DiscontGalerkinSolverCom
{

public: // functions

  /// Constructor
  explicit WallBC(const std::string& name);

  /// Destructor
  virtual ~WallBC();

  virtual void setup();

  virtual void unsetup();

  //subroutine for computation of DF matrix
  void compute_DFmatrix2D(Framework::State state, RealVector normal, RealMatrix *m_DF);
  void compute_DFmatrix3D(Framework::State state, RealVector normal, RealMatrix *m_DF);

  /// Execute processing actions
  void executeOnTrs();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  RealMatrix m_DF;

  /// handle for the InnerCells trs
  Common::SafePtr<Framework::TopologicalRegionSet> m_cells;

  /// map of LSSMatrix accumulators, one for each cell type
  Common::CFMap<CFuint,DGElemTypeData> m_mapElemData;

  /// socket for Rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

    /*                   /  invJacobi[0]  invJacobi[1]  \
    //     invJacobi =  |                                |
    //                   \  invJacobi[2]  invJacobi[3]  /
    */
  /// inverse Jacobi matrix
  //CFreal invJacobiLeft[4];

  /// determinant of Jacobi matrix
  //CFreal detJacobiLeft;

  /// pointer to temporary face in loop
  Framework::GeometricEntity * m_face;

}; // class WallBC

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_WallBC_hh

