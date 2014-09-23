#ifndef COOLFluiD_Numerics_FiniteVolume_RoeSAFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeSAFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "RoeFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux corresponding to the Euler
 * physical model 2D (in conservative variables)
 *
 * @author Andrea Lani
 *
 */
class RoeSAFlux : public RoeFlux {
public:

  /**
   * Constructor
   */
  RoeSAFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeSAFlux();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data
   */
  virtual void setup();

private: // helper function

  /**
   * Set the abs of the eigen values
   */
  virtual void setAbsEigenValues();

  /**
   * Update the etaPA parameter
   */
  void updateEtaPA(Framework::GeometricEntity *const cell, 
		   const RealVector& eValues, 
		   CFreal& etaPA);
  
private:
  
  /// array of flags for marking partition faces
  std::vector<bool> _isPartitionFace;

  /// physical model data
  RealVector _dataLeftState;

  /// physical model data
  RealVector _dataRightState;
  
  /// temporary array of eigenvalues 
  RealVector _tmpEv;
    
  /// normal vector
  RealVector _normal;

  /// ID of the entropy correction type
  CFuint _entropyFixID;
  
  /// enlarged stencil for compatibility with old COOLFluiD (wrong)
  bool _oldStencil;
  
}; // end of class RoeSAFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeSAFlux_hh
