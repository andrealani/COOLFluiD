#ifndef COOLFluiD_Numerics_MeshRigidMove_RigidMoveData_hh
#define COOLFluiD_Numerics_MeshRigidMove_RigidMoveData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "MathTools/RealVector.hh"
#include "Framework/MeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * MeshRigidMoveCom 's that compose the MeshRigidMove.
   *
   * @see MeshRigidMoveCom
   *
   * @author Thomas Wuilbaut
   */
class RigidMoveData : public Framework::MeshAdapterData {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  RigidMoveData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~RigidMoveData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up the member data
   */
   void setup();

  /**
   * Unsetup the member data
   */
  void unsetup();

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "RigidMove";
  }
  
  /**
   * Gets the Rotation Speed
   */
  CFreal getRotationSpeed()
  {
    return _rotationSpeed;
  }

  /**
   * Gets the Rotation Center
   */
  RealVector getRotationCenter()
  {
    return _rotationCenter;
  }

  /**
   * Gets the Expansion Ratio
   */
  RealVector getExpansionRatio()
  {
    return _expansionRatio;
  }

  /**
   * Gets the Rotation Center Speed
   */
  RealVector getRotationCenterSpeed()
  {
    return _rotationCenterSpeed;
  }


  /**
   * Gets the Translation Speed
   */
  RealVector getTranslationSpeed()
  {
    return _translationSpeed;
  }

private:

  /// Coordinates of center of rotation
  RealVector _rotationCenter;

  /// Speed of the center of rotation
  RealVector _rotationCenterSpeed;

  /// Translational speed
  RealVector _translationSpeed;

  /// Expansion Ratio
  RealVector _expansionRatio;

  /// Rotational Speed
  CFreal _rotationSpeed;

}; // end of class RigidMoveData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for MeshRigidMove
typedef Framework::MethodCommand<RigidMoveData> RigidMoveCom;

/// Definition of a command provider for MeshRigidMove
typedef Framework::MethodCommand<RigidMoveData>::PROVIDER RigidMoveComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshRigidMove_RigidMoveData_hh
