#include "FEM_VolumeIntegrator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

FEM_VolumeIntegrator::FEM_VolumeIntegrator() :
  VolumeIntegrator(),
  _lastPrecomputedCellID(0)
{
}

//////////////////////////////////////////////////////////////////////////////

FEM_VolumeIntegrator::~FEM_VolumeIntegrator()
{
}

//////////////////////////////////////////////////////////////////////////////

void FEM_VolumeIntegrator::precomputeCellData(Framework::GeometricEntity* geo)
{

  cf_assert(_setup);

  _intgSol = getSolutionIntegrator(geo);
  _intgGeo = getGeometryIntegrator(geo);

  cf_assert(_intgSol.isNotNull());
  cf_assert(_intgGeo.isNotNull());

  std::vector<Node*>& nodes   = *geo->getNodes();
  std::vector<State*>& states = *geo->getStates();

  _intgGeo->computeCoordinatesAtQuadraturePoints(nodes, _coord);
  const std::vector<RealVector>& mapCoord = _intgSol->getQuadraturePointsCoordinates();

  _intgGeo->computeJacobianAtQuadraturePoints(nodes, mapCoord, _jacob);

  _intgSol->computeSolutionAtQuadraturePoints(states, _solValues);

  _intgSol->computeGradSolutionShapeFAtQuadraturePoints(_jacob,mapCoord,_gradValues);

  _intgGeo->computeJacobianDetAtQuadraturePoints(nodes, _detJacobian);

  const std::valarray<CFreal>& coeff  = _intgSol->getCoeff();

  cf_assert(_solValues.size() >= coeff.size());
  cf_assert(_gradValues.size() >= coeff.size());
  cf_assert(_detJacobian.size() >= coeff.size());

//   data.solValues = &_solValues;
//   data.solShapeF = &shapeF;
//   data.coord = &_coord;
//   data.gradValues = &_gradValues;

  _lastPrecomputedCellID = geo->getID();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
