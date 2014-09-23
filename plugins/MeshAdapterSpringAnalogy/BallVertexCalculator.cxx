#include "BallVertexCalculator.hh"
#include "MathTools/MathFunctions.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

BallVertexCalculator::BallVertexCalculator() :
socket_displacements("Null"),
_normal(),
_xp(),
_ballVertex()
{

_functionMap.insert ( CFGeoShape::LINE,&BallVertexCalculator::calculateLineBallVertex);
_functionMap.insert ( CFGeoShape::TRIAG,&BallVertexCalculator::calculateTriangleBallVertex);
_functionMap.insert ( CFGeoShape::QUAD,&BallVertexCalculator::calculateQuadBallVertex);
_functionMap.insert ( CFGeoShape::TETRA,&BallVertexCalculator::calculateTetraBallVertex);
_functionMap.insert ( CFGeoShape::PYRAM,&BallVertexCalculator::calculatePyramidBallVertex);
_functionMap.insert ( CFGeoShape::PRISM,&BallVertexCalculator::calculatePrismBallVertex);
_functionMap.insert ( CFGeoShape::HEXA,&BallVertexCalculator::calculateHexaBallVertex);

_functionMap.sortKeys();

}

//////////////////////////////////////////////////////////////////////////////

BallVertexCalculator::~BallVertexCalculator()
{
}

//////////////////////////////////////////////////////////////////////////////

void BallVertexCalculator::setDataSockets(DataSocketSink<RealVector> nodalDisplacementSocket)
{
  socket_displacements = nodalDisplacementSocket;
}

//////////////////////////////////////////////////////////////////////////////

void BallVertexCalculator::setup()
{
  _normal.resize(PhysicalModelStack::getActive()->getDim());
  _xp.resize(PhysicalModelStack::getActive()->getDim());
  _ballVertex.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

RealVector BallVertexCalculator::computeBallVertex(Framework::GeometricEntity* geoEnt, CFuint iNode, CFuint jNode)
{
  _geoEnt = geoEnt;
  _iNode = iNode;
  _jNode = jNode;

  CALL_MEMBER_FN ( (*this), ( _functionMap.find ( geoEnt->getShape() ) ) ) ();

  return _ballVertex;
}

//////////////////////////////////////////////////////////////////////////////

void BallVertexCalculator::calculateTriangleBallVertex()
{

const std::vector<Node*>& nodes = *(_geoEnt->getNodes());

///Project node iNode onto the opposite face (Np)
//Determine the oposite face
CFuint kNode = 999;
if ((_iNode != 0) &&(_jNode != 0)) kNode = 0;
if ((_iNode != 1) &&(_jNode != 1)) kNode = 1;
if ((_iNode != 2) &&(_jNode != 2)) kNode = 2;
cf_assert(kNode < 3);

RealVector edgeij(PhysicalModelStack::getActive()->getDim());
RealVector edgejk(PhysicalModelStack::getActive()->getDim());
edgejk = *(nodes[kNode]) - *(nodes[_jNode]);

//Compute the normal to the oposite face
RealVector w(PhysicalModelStack::getActive()->getDim());

_normal[0] = -edgejk[1];
_normal[1] = edgejk[0];

//Normalize the normal vector
_normal.normalize();

//Determine xp
edgeij = *(nodes[_jNode]) - *(nodes[_iNode]);
w = _normal * (_normal[0] * edgeij[0] + _normal[1] * edgeij[1]);
_xp = *(nodes[_iNode]) + w;

//Interpolate displacement of virtual node Np
const CFuint nbNodes = 3;
RealVector shapeFunctions(nbNodes);
shapeFunctions = _geoEnt->computeShapeFunctionAtCoord(_xp);
DataHandle<RealVector> displacements = socket_displacements.getDataHandle();

_ballVertex = 0.;
for(CFuint iNode=0; iNode < nbNodes; iNode++)
{
  CFuint localID = nodes[iNode]->getLocalID();
  _ballVertex += displacements[localID] * shapeFunctions[iNode];

}

// if(_ballVertex.sum()>0){
//
// std::cout << "DisplacementBallvertex: " << _ballVertex << std::endl;
// std::cout << "Displacement[0]: " << _displacements[nodes[0]->getLocalID()] << std::endl;
// std::cout << "Displacement[1]: " << _displacements[nodes[1]->getLocalID()] << std::endl;
// std::cout << "Displacement[2]: " << _displacements[nodes[2]->getLocalID()] << std::endl;
// std::cout << "Node ["<< _iNode<<"]: " << *(nodes[_iNode]) << std::endl;
// std::cout << "Node ["<< _jNode<<"]: " << *(nodes[_jNode]) << std::endl;
// std::cout << "Node ["<< kNode<<"]: " << *(nodes[kNode]) << std::endl;
// std::cout << "Projected point" << _xp << std::endl;
// std::cout << "ShapeFunctions: " << shapeFunctions << std::endl;
// std::cout << "----------------------------------------------------------------" << std::endl;
// }

}

//////////////////////////////////////////////////////////////////////////////

void BallVertexCalculator::calculateTetraBallVertex()
{
 throw Common::NotImplementedException (FromHere(),"BallVertexCalculator::calculateTetraBallVertex()");

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

