#ifndef COOLFluiD_UFEM_functors_GradSquared_hh
#define COOLFluiD_UFEM_functors_GradSquared_hh

#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"

namespace COOLFluiD {

namespace UFEM {

namespace functors {

using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

/// Evaluates the square of the gradient of the shape functions for the given mapped coordinates
struct GradSquared
{
  typedef COOLFluiD::RealMatrix value_t;
  /// Cell associated with this functor
  const GeometricEntity& cell;

  /// Constants for the coordinate interpolation
  CFreal ax, bx, cx, dx, ay, by, cy, dy;

  /// Construct using a GeometricEntity
  /**
   * @param GE GeometricEntity overwhich the evaluation is done
   */
  GradSquared(const GeometricEntity& GE) : cell(GE) {
    const std::vector<Node*>& nodes = cell.getNodes();
    ax =  nodes[0]->at(XX) + nodes[1]->at(XX) + nodes[2]->at(XX) + nodes[3]->at(XX);
    bx = -nodes[0]->at(XX) + nodes[1]->at(XX) + nodes[2]->at(XX) - nodes[3]->at(XX);
    cx = -nodes[0]->at(XX) - nodes[1]->at(XX) + nodes[2]->at(XX) + nodes[3]->at(XX);
    dx =  nodes[0]->at(XX) - nodes[1]->at(XX) + nodes[2]->at(XX) - nodes[3]->at(XX);
    ay =  nodes[0]->at(YY) + nodes[1]->at(YY) + nodes[2]->at(YY) + nodes[3]->at(YY);
    by = -nodes[0]->at(YY) + nodes[1]->at(YY) + nodes[2]->at(YY) - nodes[3]->at(YY);
    cy = -nodes[0]->at(YY) - nodes[1]->at(YY) + nodes[2]->at(YY) + nodes[3]->at(YY);
    dy =  nodes[0]->at(YY) - nodes[1]->at(YY) + nodes[2]->at(YY) - nodes[3]->at(YY);
  }

  /// Returns the total number of states for the cell
  CFuint size()
  {
    return cell.nbStates()*cell.nbStates();
  }

  /// Evaluate function for the given real coordinates
  const COOLFluiD::RealMatrix operator()(const RealVector& mapped_coord)
  {
    std::vector<RealVector> mapped_coords;
    mapped_coords.push_back(mapped_coord);
    GeometricEntity& non_const_ge = const_cast<GeometricEntity&>(cell);
    COOLFluiD::RealMatrix B = non_const_ge.computeSolutionShapeFunctionGradients(mapped_coords).front();
    COOLFluiD::RealMatrix B_trans(B.nbCols(), B.nbRows());
    B.transpose(B_trans);
    return B*B_trans;
  }

  /// Evalate the function in the given mapped coordinates xi and eta
  const COOLFluiD::RealMatrix operator()(const CFreal xi, const CFreal eta)
  {
    RealMatrix grad_mapped(4, 2); // Stores the gradient in mapped coordinates
    grad_mapped(0, XX) = -1 + eta;
    grad_mapped(1, XX) =  1 - eta;
    grad_mapped(2, XX) =  1 + eta;
    grad_mapped(3, XX) = -1 - eta;
    grad_mapped(0, YY) = -1 + xi;
    grad_mapped(1, YY) = -1 - xi;
    grad_mapped(2, YY) =  1 + xi;
    grad_mapped(3, YY) =  1 - xi;
    RealMatrix jac_adj(2, 2); // Adjoint of the jacobian
    jac_adj(0, XX) =  cy + dy*xi;
    jac_adj(1, XX) = -by - dy*eta;
    jac_adj(0, YY) = -cx - dx*xi;
    jac_adj(1, YY) =  bx + dx*eta;

    const CFreal jac_det = (bx*dy - by*dx)*xi + (cy*dx - cx*dy)*eta + bx*cy - by*cx;
    COOLFluiD::RealMatrix B(4, 2);
    B = grad_mapped * jac_adj;
    COOLFluiD::RealMatrix B_trans(B.nbCols(), B.nbRows());
    B.transpose(B_trans);
    return B*B_trans/(16.*jac_det);
  }

  /// Evaluate the function in mapped coordinates
  /**
   * @param node A node storing the mapped coordinates.
   */
  const COOLFluiD::RealMatrix operator()(Node* node)
  {
    return operator()(*node->getData());
  }
}; // struct GradSquared

} // namespace functors

} // namespace UFEM

} // namespace COOLFluiD

#endif // COOLFluiD_UFEM_functors_GradSquared_hh$
