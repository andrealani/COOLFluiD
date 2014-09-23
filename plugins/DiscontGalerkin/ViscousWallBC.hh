#ifndef COOLFluiD_DiscontGalerkin_ViscousWallBC_hh
#define COOLFluiD_DiscontGalerkin_ViscousWallBC_hh

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
 * This is a standard command to assemble the system using
 * discontinuous Galerkin solver by adding the contributions integrals
 * in the faces.
 *
 * @author Martin Holik
 */
class ViscousWallBC : public ViscousBaseSolve
{

public: // functions

  /// Constructor
  explicit ViscousWallBC(const std::string& name);

  /// Destructor
  virtual ~ViscousWallBC();

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

  /// matrix of linearization of viscous terms K[s][k](x,x)
  std::vector < std::vector < RealMatrix > > m_kMatrix;

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

  // constant for viscous flow
  CFreal m_sigma;
  CFreal m_theta;

}; // class ViscousWallBC

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_ViscousWallBC_hh

