#ifndef COOLFluiD_Numerics_FiniteVolume_DiamondVolume3DDerivative_hh
#define COOLFluiD_Numerics_FiniteVolume_DiamondVolume3DDerivative_hh

//////////////////////////////////////////////////////////////////////////////

#include "DerivativeComputer.hh"
#include "Framework/NormalsCalculator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes derivatives using a diamond volume in 3D
 *
 * @author Andrea Lani
 *
 */
class DiamondVolume3DDerivative : public DerivativeComputer {
public:

  /**
   * Constructor
   */
  DiamondVolume3DDerivative(const std::string& name);

  /**
   * Default destructor
   */
  ~DiamondVolume3DDerivative();

  /**
   * Set up the member data
   */
  virtual void setup();

  /*
   * Compute the gradients
   */
  void computeGradients(const RealMatrix& values,
			std::vector<RealVector*>& gradients);

  /**
   * Compute the control volume around the current face
   */
  void computeControlVolume(std::vector<RealVector*>& states,
			    Framework::GeometricEntity *const geo);
  
  /**
   * Compute the average values corresponding to the given values
   */
  void computeAverageValues(Framework::GeometricEntity *const geo,
			    const std::vector<RealVector*>& values,
			    RealVector& avValues);
  
  /**
   * Get the maximum number of vertices in the control volume
   */
  CFuint getMaxNbVerticesInControlVolume() const
  {
    return 6;
  }

  /**
   * Get the current number of vertices in the control volume
   */
  CFuint getNbVerticesInControlVolume(Framework::GeometricEntity *const geo) const;
  
  /**
   * Get the jacobian of the gradients
   */
  Common::SafePtr<std::vector<RealVector> > getGradientsJacob();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
    
private: // methods

  /**
   * Compute 2 tetras volume
   */
  void compute2TetraVolume(Framework::GeometricEntity *const geo);

  /**
   * Compute 2 pyramids volume
   */
  void compute2PyramidVolume(Framework::GeometricEntity *const geo);
  
  /**
   * Compute 2 pyramids normals
   */
  void compute2PyramidNormals();
  
  /**
   * Compute 2 tetras gradient
   */
  void compute2TetraGradient(const RealMatrix& values,
			     std::vector<RealVector*>& gradients);

  /**
   * Compute 2 pyramids gradient
   */
  void compute2PyramidGradient(const RealMatrix& values,
			       std::vector<RealVector*>& gradients);
  
  /**
   * Planarize the face
   * @post nodes will be repositioned
   */
  void planarizeFaceStart(Framework::GeometricEntity& geo); 
  
  /**
   * Restore settings as before planarization
   */
  void planarizeFaceEnd(Framework::GeometricEntity& geo)
  {
    using namespace COOLFluiD::Framework;
    
    if (_nonPlanarFaceNodeID[geo.getID()] != -1) {
      *geo.getNode(_nonPlanarFaceNodeID[geo.getID()]) = _ncoord;
      geo.getNeighborGeo(0)->getState(0)->getCoordinates() = _lcoord;
      if (!geo.getState(1)->isGhost()) {
	geo.getNeighborGeo(1)->getState(0)->getCoordinates() = _rcoord;
      }
    }
  }
  
  /**
   * Correct non planar face
   */
  void correctNonPlanarFace(Framework::GeometricEntity& geo, CFint faceNodeID)
  { 
    using namespace COOLFluiD::Framework;
    Node& node = *geo.getNode(faceNodeID);
    _ncoord = node;
    node = _xproj;
    _nonPlanarFaceNodeID[geo.getID()] = faceNodeID;
  }
  
  /**
   * Recompute the cell centers
   */
  void recomputeCellCenter(Framework::GeometricEntity *const cell, RealVector& bkpcc)
  {
    const CFuint nbNodes = cell->nbNodes();
    const CFreal ovNbNodes = 1./static_cast<CFreal>(nbNodes);
    RealVector& ccenter = cell->getState(0)->getCoordinates();
    bkpcc = ccenter;
    ccenter = 0.;
    for (CFuint i = 0; i < nbNodes; ++i) {
      ccenter += ovNbNodes*(*cell->getNode(i));
      //std::cout.precision(14); std::cout.setf(std::ios::scientific,std::ios::floatfield); std::cout << *cell->getNode(i) << std::endl;;
    }
    //std::cout << std::endl;
  }  
  
private: // data

  /// one third
  const CFreal _third;

  /// storage of states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// storage of nodes
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of cell flags
  Framework::DataSocketSink<bool> socket_cellFlag;
  
  /// normals calculator
  Framework::NormalsCalculator _normalCalculator;

  /// normal 031 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n031;

  /// normal 132 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n132;

  /// normal 023 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n023;

  /// normal 041 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n041;

  /// normal 142 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n142;

  /// normal 024 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n024;

  /// mid node coordinates
  RealVector _midNode;
  
  /// temporary normal
  RealVector _tmpNormal;
  
  /// array
  RealVector _v1;

  /// array
  RealVector _v2;

  /// array
  RealVector _v3;

  /// array
  RealVector _v4;

  /// array
  RealVector _xproj;
  
  /// array
  RealVector _ncoord;
  
  /// array
  RealVector _lcoord;
  
  /// array
  RealVector _rcoord;
  
  /// tetra coordinates
  RealMatrix _tetraCoord;

  /// pyramid coordinates
  RealMatrix _pyramCoord;

  /// pyramid normals
  RealMatrix _pyramNormals1;

  /// pyramid normals
  RealMatrix _pyramNormals2;
  
  // volume centroid
  RealVector _cVolume;
  
  // face centroid
  RealVector _cFace;

  // vector face normal - face centroid
  RealVector _ncF;
  
  // vector volume centroid - face centroid
  RealVector _cVcF;
  
  // vector vertex pyramid - face centroid
  RealVector _n4cF;
  
  // vector vertex pyramid - face centroid
  RealVector _n5cF;
  
  /// control volumes for each face
  RealVector _cvolumes;
  
  /// flag storage telling if normal must be changed of sign
  std::vector<bool> _changeNormalSign;
  
  /// ID of the non planar node in each face
  std::vector<CFint> _nonPlanarFaceNodeID;
  
}; // end of class DiamondVolume3DDerivative

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DiamondVolume3DDerivative_hh
