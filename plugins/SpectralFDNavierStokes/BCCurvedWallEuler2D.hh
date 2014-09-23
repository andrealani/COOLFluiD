#ifndef COOLFluiD_Numerics_SpectralFD_BCCurvedWallEuler2D_hh
#define COOLFluiD_Numerics_SpectralFD_BCCurvedWallEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/VectorialFunction.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a curved wall boundary condition, using the Krivodonova and Berger treatment.
 * It is meant to be used in combination with P1 elements.
 * @see L. Krivodonova and M. Berger, "High-order accurate implementation of solid wall
 *      boundary conditions in curved geometries", J. Comput. Phys. 211 (2006) 492-512
 *
 * @author Kris Van den Abeele
 */
class BCCurvedWallEuler2D : public BCStateComputer {

public:  // methods

  /**
 * Defines the Config Option's of this class
 * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCCurvedWallEuler2D(const std::string& name);

  /// Destructor
  ~BCCurvedWallEuler2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCCurvedWallEuler2D";
  }

  /// Set up private data and data
  void setup();

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                             std::vector< std::vector< RealVector* > >& ghostGrads,
                             const std::vector< RealVector >& normals,
                             const std::vector< RealVector >& coords);

protected: // methods

  /**
   * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );

private: // methods

  /**
   * Computes the normals in the face flux points from mesh data
   */
  void computeNormalsInFaceFluxPointsFromUserFunction();

  /**
   * Computes the normals in the face flux points from mesh data
   */
  void computeNormalsInFaceFluxPointsFromMeshData();

  /**
   * Computes the normals in the face flux points from the domain model
   */
  void computeNormalsInFaceFluxPointsFromDomainModel();

  /**
   * Sets the derivative of the shape functions at the given coordinate
   */
  void setShapeFunctionDerivatives(const CFreal ksi, RealVector& shapeFuncDerivs);

protected: // data

  /// handle to the Node's storage
  Framework::DataHandle< Framework::Node*, Framework::GLOBAL> m_nodes;

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the variables
  std::vector<std::string> m_vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// flux point local coordinates on the faces
  Common::SafePtr< std::vector< RealVector > > m_flxPntFaceLocalCoords;

  /// map of boundary face global ID's to normals in the face flux points
  std::map< CFuint , std::vector< RealVector > > m_faceFlxPntNormals;

  /// boolean telling whether to use the domain model for the computation of the normals
  bool m_useDomainModel;

}; // class BCCurvedWallEuler2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BCCurvedWallEuler2D_hh
