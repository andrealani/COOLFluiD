#ifndef COOLFluiD_DiscontGalerkin_SuperOutletBC_hh
#define COOLFluiD_DiscontGalerkin_SuperOutletBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Framework/DataSocketSink.hh"
#include "DiscontGalerkin/DGElemTypeData.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"
#include "DiscontGalerkin/StdBaseSolve.hh"
#include "Framework/VectorialFunction.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to assemble the system using
 * discontinuous Galerkin solver by adding the contributions integrals
 * in the boundary faces.
 *
 * @author Martin Holik
 */
class SuperOutletBC : public StdBaseSolve
{

public: // functions

  /// Constructor
  explicit SuperOutletBC(const std::string& name);

  /// Destructor
  virtual ~SuperOutletBC();

  virtual void setup();

  virtual void unsetup();

  /// Execute processing actions
  void executeOnTrs();

 /**
  * Returns the DataSocket's that this command needs as sinks
  * @return a vector of SafePtr with the DataSockets
  */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

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

  /// pointer to temporary face in loop
  Framework::GeometricEntity * m_face;

  CFreal kappa1;

}; // class SuperOutletBC

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_SuperOutletBC_hh

