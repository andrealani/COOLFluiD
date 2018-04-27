#ifndef COOLFluiD_Numerics_FiniteVolume_NeumannBC_hh
#define COOLFluiD_Numerics_FiniteVolume_NeumannBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NormalsCalculator.hh"
#include "Framework/VolumeCalculator.hh"
#include "FiniteVolume/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Neumann boundary condition
   *
   * @author Andrea Lani
   *
   */
class NeumannBC : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NeumannBC(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NeumannBC();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();  

protected:

  /// Compute the 2D control volume around the current face
  void computeControlVolume2D(std::vector<RealVector*>& states, 
			      Framework::GeometricEntity *const geo);

  /// Compute the 3D control volume made by 2 tetrahedra
  void compute2TetraVolume(std::vector<RealVector*>& states, 
			   Framework::GeometricEntity *const geo);
  
  /// Compute the 3D control volume made by 2 pyramids
  void compute2PyramVolume(std::vector<RealVector*>& states, 
			   Framework::GeometricEntity *const geo);
  
  /// Compute the ghost state with only radial gradient activated in 2D 
  void computeGhostWithRadialGradient2D(Framework::GeometricEntity *const face);
  
  /// Compute the ghost state with only radial gradient activated in 3D 
  void computeGhostWithRadialGradient3DTetra(Framework::GeometricEntity *const face);
  
  /// Compute the ghost state with only radial gradient activated in 3D 
  void computeGhostWithRadialGradient3DPyram(Framework::GeometricEntity *const face);
  
private: // data

  /// storage of the nodal states
  Framework::DataSocketSink<RealVector> socket_nstates;

  /// volume calculator
  Framework::VolumeCalculator _volumeCalculator;
  
  /// normals calculator
  Framework::NormalsCalculator _normalCalculator;
  
  /// normal 01 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n01;

  /// normal 12 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n12;
  
  /// normal 23 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n23;
  
  /// normal 30 in the quadrilateral control volume for
  /// gradients evaluation
  RealVector _n30;
  
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
  
  // array of states involved in the flux computation
  std::vector<RealVector*> _states;
  
  /// quad control volume nodes
  std::vector<Framework::Node*> _nodes;
    
  /// null the gradients in theta and phi
  bool _onlyRadialGradient;
  
}; // end of class NeumannBC

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NeumannBC_hh
