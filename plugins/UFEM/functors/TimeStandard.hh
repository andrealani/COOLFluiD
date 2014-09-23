#ifndef COOLFluiD_UFEM_functors_TimeStandard_hh
#define COOLFluiD_UFEM_functors_TimeStandard_hh

#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"

namespace COOLFluiD {

namespace UFEM {

namespace functors {

using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

/// Evaluates the time-dependent matrix contribution for standard FEM
struct TimeStandard
{
  typedef COOLFluiD::RealMatrix value_t;
  const GeometricEntity& cell;

  /// Construct using a GeometricEntity
  /**
   * @param GE GeometricEntity overwhich the evaluation is done
   */
  TimeStandard(const GeometricEntity& GE) : cell(GE) {
    CFout << "------------ cell: " << GE.getID() << "---------------\n";
  }

  /// Returns the total number of states for the cell
  CFuint size()
  {
    return cell.nbStates()*cell.nbStates();
  }

  const COOLFluiD::RealMatrix operator()(const RealVector& coord)
  {
    // todo: const-correctness
    GeometricEntity& cast_cell = const_cast<GeometricEntity&>(cell);
    COOLFluiD::RealVector shape_funcs = cast_cell.computeShapeFunctionAtMappedCoord(coord);
    const CFuint n = shape_funcs.size();
    COOLFluiD::RealMatrix N(n, 1);
    for(CFuint i = 0; i != n; ++i) {
      N(i, 0) = shape_funcs[i];
    }
    COOLFluiD::RealMatrix N_trans(1, n);
    N.transpose(N_trans);
    COOLFluiD::RealMatrix result = N * N_trans;
    CFout << "time result: " << result << "\n";
    std::vector<RealVector> mapped_coords;
    mapped_coords.push_back(coord);
    return N * N_trans * cell.computeGeometricShapeFunctionJacobianDeterminant(mapped_coords)[0];
  }

  /// Evalate the function in the given mapped coordinates xi and eta
  const COOLFluiD::RealMatrix operator()(const CFreal xi, const CFreal eta)
  {
    RealVector mapped_coord(2);
    mapped_coord[KSI] = xi;
    mapped_coord[ETA] = eta;
    return operator()(mapped_coord);
  }

  /// Evaluate the function in mapped coordinates
  /**
   * @param node A node storing the mapped coordinates.
   */
  const COOLFluiD::RealMatrix operator()(Node* node)
  {
    return operator()(*node->getData());
  }
}; // struct TimeStandard

} // namespace functors

} // namespace UFEM

} // namespace COOLFluiD

#endif // COOLFluiD_UFEM_functors_TimeStandard_hh$
