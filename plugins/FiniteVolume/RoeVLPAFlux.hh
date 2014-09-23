#ifndef COOLFluiD_Numerics_FiniteVolume_RoeVLPAFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeVLPAFlux_hh

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
class RoeVLPAFlux : public RoeFlux {
public:

  /**
   * Constructor
   */
  RoeVLPAFlux(const std::string& name);

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Default destructor
   */
  virtual ~RoeVLPAFlux();

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
   * Compute the etaVL parameter
   */
  void computeEtaVL(CFreal& etaVL) const
  {
    etaVL = 0.0;	
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i < nbEqs; ++i) {
      etaVL = std::max(etaVL, _rightEvalues[i] - _leftEvalues[i]);
      // etaVL = std::max(etaVL,std::abs(_rightEvalues[i] - _leftEvalues[i]));
    }	
    //etaVL *= 0.5;	
    
    //  const CFuint vlSize = _vlFixIDs.size();
    // for (CFuint i = 0; i < vlSize; ++i) {
    //   etaVL = std::max(etaVL, _rightEvalues[_vlFixIDs[i]] - _leftEvalues[_vlFixIDs[i]]);
    // }
  }
  
  /**
   * Update the etaPA parameter
   */
  void updateEtaPA(Framework::GeometricEntity *const cell, 
		   const RealVector& eValues, 
		   Framework::State *const otherState,
		   CFreal& etaPA);
  
private:

  /// physical model data
  RealVector _dataLeftState;

  /// physical model data
  RealVector _dataRightState;
  
  /// temporary array of eigenvalues 
  RealVector _tmpEv;
    
  /// normal vector
  RealVector _normal;

  /// use Pandolfi - D'Ambrosio's fix
  bool _useFixPA;
  
  /// array telling the IDs of the eigenvalues on which VL fix must be applied
  std::vector<CFuint> _vlFixIDs;
  
  /// array of flags to tell if VL fix must be applied
  std::vector<bool> _flagVL;

}; // end of class RoeVLPAFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeVLPAFlux_hh
