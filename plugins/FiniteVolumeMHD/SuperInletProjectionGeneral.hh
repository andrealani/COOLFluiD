#ifndef COOLFluiD_Numerics_FiniteVolume_SuperInletProjectionGeneral_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperInletProjectionGeneral_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class MapGeoToTrsAndIdx;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic inlet command in 2D
   * for projection scheme
   *
   * @author Radka Keslerova
   *
   */
class SuperInletProjectionGeneral : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SuperInletProjectionGeneral(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SuperInletProjectionGeneral();

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
  
private: //data

  /// pointer to the mapping face - TRS
  Common::SafePtr<Framework::MapGeoToTrsAndIdx> m_mapGeoToTrs;
  
  /// phi value that is to be fixed
  bool   m_BfromFile;
  CFint _inletCoronalBC;
  CFint _Phi_divB_zero;
  CFint _Phi_divB_extrapolated;
  CFint _hydrodynamic_limit;
  CFint _DifferentialRotation;
  CFint _pressure_fixed;
  CFint _pressure_Neumann;
  CFreal _pBC;
  CFreal _rhoBC;
  CFreal _VrBC;  
  CFint _rotation;
  CFreal _betamin;
  CFreal _vAmax;
  bool _nonhomogeneous; 

  /// array specifying the IDs for which a special treatment has to be applied
  std::vector<CFuint> _projectionIDs;
  
  /// IDs of the variables from which values are read by file
  std::vector<CFuint> m_varIDs;
  
}; // end of class SuperInletProjectionGeneral

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperInletProjectionGeneral_hh
