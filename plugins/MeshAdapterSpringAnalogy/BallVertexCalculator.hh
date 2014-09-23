#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_BallVertexCalculator_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_BallVertexCalculator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/Storage.hh"
#include "Common/CFMap.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class represents a calculator of the BallVertex in elements
 *
 * @author Thomas Wuilbaut
 *
 */

class BallVertexCalculator {

/// pointer to member function
typedef void (BallVertexCalculator::*Computer)();

public:

  /**
   * Default constructor without arguments.
   */
  BallVertexCalculator();

  /**
   * Default destructor.
   */
  ~BallVertexCalculator();

  void setup();

  /**
   * Calculate the BallVertex of a cell
   * @param coord coordinates of the nodes
   */
  RealVector computeBallVertex(Framework::GeometricEntity* geoEnt, CFuint iNode, CFuint jNode);

private: // functions

  /**
   * Compute BallVertex in the case of a Line
   */
  void calculateLineBallVertex()
  {
    throw Common::NotImplementedException (FromHere(),"BallVertexCalculator::calculateLineBallVertex()");
  }

  /**
   * Compute BallVertex in the case of a Line
   */
  void calculateTriangleBallVertex();

  /**
   * Compute BallVertex in the case of a Line
   */
  void calculateQuadBallVertex()
  {
    throw Common::NotImplementedException (FromHere(),"BallVertexCalculator::calculateQuadBallVertex()");
  }

  /**
   * Compute BallVertex in the case of a Line
   */
  void calculateTetraBallVertex();

  /**
   * Compute BallVertex in the case of a Line
   */
  void calculatePrismBallVertex()
  {
    throw Common::NotImplementedException (FromHere(),"BallVertexCalculator::calculatePrismBallVertex()");
  }

  /**
   * Compute BallVertex in the case of a Line
   */
  void calculatePyramidBallVertex()
  {
    throw Common::NotImplementedException (FromHere(),"BallVertexCalculator::calculatePyramidBallVertex()");
  }

  /**
   * Set the data sockets to be used by the BallVertex Calculator
   */
  void setDataSockets(Framework::DataSocketSink<RealVector> nodalDisplacementSocket);

  /**
   * Compute BallVertex in the case of a Line
   */
  void calculateHexaBallVertex()
  {
    throw Common::NotImplementedException (FromHere(),"BallVertexCalculator::calculateHexaBallVertex()");
  }

private: // data

  /// socket for nodal displacement
  Framework::DataSocketSink<RealVector > socket_displacements;

  Common::CFMap<CFGeoShape::Type,Computer> _functionMap;

  //idx of the iNode
  CFuint _iNode;

  // idx of the j node
  CFuint _jNode;

  // the normal
  RealVector _normal;

  // the coordinate of the projected point
  RealVector _xp;

  // the resulting ballVertex spring
  RealVector _ballVertex;

  // the element
  Framework::GeometricEntity* _geoEnt;

}; // end of class BallVertexCalculator

//////////////////////////////////////////////////////////////////////////////
    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshAdapterSpringAnalogy_BallVertexCalculator_hh
