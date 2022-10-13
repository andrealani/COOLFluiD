#ifndef COOLFluiD_Numerics_FiniteVolume_SuperInletProjection_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperInletProjection_hh

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
class SuperInletProjection : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SuperInletProjection(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SuperInletProjection();

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
  
  /// map each wall faceID to the corresponding wall temperature
  Common::CFMap<Framework::TopologicalRegionSet*, RealVector*> m_mapTrs2Twall;

  /// phi value that is to be fixed
  bool   m_BfromFile;
  CFreal _refPhi;
  CFint _inletCoronalBC;
  CFint _Phi_divB_zero;
  CFint _Phi_divB_extrapolated;
  CFint _JensVelocityBC;
  CFint _BarbarasVelocityBC;
  CFint _hydrodynamic_limit;
  CFint _DanasVelocityBC;
  CFint _DifferentialRotation;
  CFint _JensBfieldBC;
  CFint _DanasBfieldBC;
  CFint _JonLinkersBfieldSuggestion;
  CFint _pressure_fixed;
  CFint _pressure_Neumann;
  CFint _JensRhoIni;
  CFint _JensPIni;
  CFreal _pBC;
  CFreal _rhoBC;
  CFreal _VrBC;  
  CFint _rotation;
  /// array specifying the IDs for which a special treatment has to be applied
  std::vector<CFuint> _projectionIDs;
  
  /// IDs of the variables from which values are read by file
  std::vector<CFuint> m_varIDs;
  
}; // end of class SuperInletProjection

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperInletProjection_hh
