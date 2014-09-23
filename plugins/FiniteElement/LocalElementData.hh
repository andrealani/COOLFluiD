#ifndef COOLFluiD_Numerics_FiniteElement_LocalElementData_hh
#define COOLFluiD_Numerics_FiniteElement_LocalElementData_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElement/FElemTypeData.hh"
#include "Framework/TopologicalRegionSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/// Data of the current element being processed to be shared between all
/// commands and strategies of the FiniteElement method
class LocalElementData
{
  public: // data

  //Cell Related Data

  /// pointer to the current Trs being processed
  Common::SafePtr<Framework::TopologicalRegionSet> trs;

  /// pointer to the current cell being processed
  Framework::GeometricEntity* cell;

  /// number of states in this element
  CFuint nbStates;

  /// cache the number of equations
  CFuint nbEqs;

  /// BlockAccumulator for the current element
  Framework::BlockAccumulator* blockacc;

  /// Local stiffness matrix
  RealMatrix* stiff_mat;

  /// Local load vector
  RealVector* load_vec;

  /// Nodal residuals
  std::vector<RealVector>* residual;

  //State Related Data

  /// state relative to the rows of the matrix
  CFuint iState;

  /// state relative to the rows of the matrix
  CFuint jState;

  //Integration Related Data
  ///Current Quadrature Point
  CFuint quadPointID;

  ///Solutions Value at Quadrature Points
  std::vector<Framework::State*>* solValues;

  ///Solution Shape Function at Quadrature Points
  const std::vector<RealVector>* solShapeF;

  /// Coordinates of the Quadrature Points
  std::vector<Framework::Node*>* coord;

  /// storage for temporary interpolated gradient values
  /// size of vector is maximun size of quadrature points
  /// each matris is sized nb shape functions * nb of dimensions
  std::vector<RealMatrix>* gradValues;

}; // end LocalElementData

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_LocalElementData_hh




