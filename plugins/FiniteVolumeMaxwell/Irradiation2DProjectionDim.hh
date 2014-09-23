#ifndef COOLFluiD_Numerics_FiniteVolume_Irradiation2DProjectionDim_hh
#define COOLFluiD_Numerics_FiniteVolume_Irradiation2DProjectionDim_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace Maxwell {
      class Maxwell2DProjectionAdimVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a irradiation inlet command
   *
   * @author Alejandro Alvarez
   *
   */
class Irradiation2DProjectionDim : public SuperInlet {
public:

  /**
   * Constructor
   */
  Irradiation2DProjectionDim(const std::string& name);

  /**
   * Default destructor
   */
  ~Irradiation2DProjectionDim();

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
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

protected: // data

  /// physical model var set
  Common::SafePtr<Physics::Maxwell::Maxwell2DProjectionAdimVarSet> _varSet;
  
  /// storage for the temporary boundary point coordinates
  RealVector _variables;

}; // end of class Irradiation2DProjectionDim

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Irradiation2DProjectionDim_hh
