#ifndef COOLFluiD_Numerics_FiniteVolume_RoeSAFluxGhost_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeSAFluxGhost_hh

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
class RoeSAFluxGhost : public RoeFlux {
public:

  /**
   * Constructor
   */
  RoeSAFluxGhost(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeSAFluxGhost();

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
  
  /**
   * Update the etaPA parameter
   */
  void updateEtaPAGhost(const RealVector& lEvalues, 
			const RealVector& rEvalues, 
			CFreal& etaPA);
  
private:
  
  // builder of cells from TRS 
  Framework::CellTrsGeoBuilder m_cellBuilder;
  
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
  
}; // end of class RoeSAFluxGhost

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeSAFluxGhost_hh
