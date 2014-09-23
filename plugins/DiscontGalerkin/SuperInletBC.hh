#ifndef COOLFluiD_DiscontGalerkin_SuperInletBC_hh
#define COOLFluiD_DiscontGalerkin_SuperInletBC_hh

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
class SuperInletBC : public StdBaseSolve
{

public: // functions

  /// Constructor
  explicit SuperInletBC(const std::string& name);

  /// Destructor
  virtual ~SuperInletBC();

  virtual void setup();

  virtual void unsetup();

  /// Execute processing actions
  void executeOnTrs();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

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

  /// determinant of Jacobi matrix
  //CFreal detJacobiLeft;

  /// primitive shape function
  //std::vector<CFreal> *m_shapeFaceFunction[3][3][3];//[function][derivation][face][point of quadrature]

  ///matrix of multiplication of shape function and its derivation on face of elements i and j
  //std::vector<CFreal> *m_shapeFaceMatrix[3][3][3][3][3]; //phi_i *phi_j   [i][j][derivation of i][derivation of j][face on i][point of quadrature]

  /// pointer to temporary face in loop
  Framework::GeometricEntity * m_face;

  CFreal kappa1;

protected:

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

}; // class SuperInletBC

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_SuperInletBC_hh

