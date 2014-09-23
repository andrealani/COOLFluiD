#ifndef COOLFluiD_Numerics_SubSytemCoupler_NewtonMeshMatcherWrite_hh
#define COOLFluiD_Numerics_SubSytemCoupler_NewtonMeshMatcherWrite_hh

//////////////////////////////////////////////////////////////////////////////

#include "StdMeshMatcherWrite.hh"
#include "Framework/NumericalJacobian.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to set the
   * match between meshes
   *
   * @author Thomas Wuilbaut
   *
   */

class NewtonMeshMatcherWrite : public StdMeshMatcherWrite {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit NewtonMeshMatcherWrite(const std::string& name);

  /**
   * Destructor.
   */
  ~NewtonMeshMatcherWrite();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

protected: // data

  /**
   * Executes the command.
   */
  virtual void executeWrite(const CFuint iProc);

  /**
   * Preselect the faces on which to apply the projection algorithm
   */
  virtual void facePreSelection();

  /**
   * Pairs the projected point with the closest face
   */
  virtual void nodeToElementPairing(SubSysCouplerData::GeoEntityIdx& matchingFace);

  /**
   * Projects the point on the closest face
   */
  virtual void computeProjection();

  /**
   * Newton Loop: Setup and Initialize the vectors
   */
  virtual void setupNewton();

  /**
   * Newton Loop: Computes the residual and the jacobian
   */
  virtual void takeStep();

  /**
   * Newton Loop: Updates the Solution
   */
  virtual void updateSolution();

  /**
   * Modified Newton Loop: Updates the Solution
   */
  virtual void newUpdateSolution();

  /**
   * Newton Loop: Computes the Residual
   */
  virtual CFreal computeResidual();

private:

void unshiftStatesCoord();
void shiftStatesCoord();

protected: // data

  /// Numerical jacobian calculator
  std::auto_ptr<Framework::NumericalJacobian> _numericalJacob;

  /// Pointer to the current face
  Framework::GeometricEntity * _currentFace;

  /// Mapped Coordinates
  RealVector _mappedCoord;

  /// Mapped Coordinates for the projection
  RealVector _mappedCoordMin;

  /// Coordinates of the point to be projected
  RealVector _coord;

  /// Temporary Coordinates Vector
  RealVector _tempCoord;

  /// Reference values for the numerical jacobian
  RealVector _refValues;

  /// Coordinates of the centroid of the faces
  RealVector _centroid;

  /// Residual
  RealVector _residual;

  /// Perturbed Residual
  RealVector _otherResidual;

  /// Jacobian
  RealVector _jacobian;

  /// Selection of faces with the closest centroid
  std::vector<SubSysCouplerData::GeoEntityIdx> _faceSelection;

  /// iVar to be modified for the step
  CFuint _currentIVar;

  /// maximum accuracy for finding the projection
  CFreal _maxAccuracy;

  /// maximum number of steps for finding the projection
  CFuint _maxNewtonSteps;

  /// number of preselected (based on centroid distance) faces
  CFuint _nbSelectedFaces;

  RealVector _bestCoord;
  CFreal _bestDistance;

  bool _shiftStates;

}; // class NewtonMeshMatcherWrite

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_NewtonMeshMatcherWrite_hh

