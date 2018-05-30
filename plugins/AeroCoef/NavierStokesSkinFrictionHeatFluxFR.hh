#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxFR_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxFR_hh

//////////////////////////////////////////////////////////////////////////////

#include "AeroCoef/AeroForcesFR.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }

  namespace FluxReconstructionMethod {
    class FluxReconstructionElementData;
  }
  
  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the skin friction and the heat flux for NavierStokes
 * simulations with FluxReconstruction
 *
 * @author Ray Vandenhoeck
 *
 */
class NavierStokesSkinFrictionHeatFluxFR : public AeroForcesFR {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesSkinFrictionHeatFluxFR(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesSkinFrictionHeatFluxFR();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();

protected:
  
  /**
   * Open the Output File and Write the header
   */
  virtual void prepareOutputFileWall();
  
  /**
   * Update the Output file with the wall values
   */
  virtual void updateOutputFileWall();
  
  /**
   * Compute the required values
   */
  virtual void computeWall();
  
  /**
   * Compute the required values
   */
  virtual void computeExtraValues()
  {
  }
  
  /**
   * Compute the value of Y+
   */
  void computeYplus();
  
  /// Update the data to write
  virtual void updateWriteData(CFuint flxIdx);
  
  /**
   * Compute the value of tau at the wall
   */
  void computeTauWall(CFuint flxIdx);
  
  /**
   * Compute dimensional pressure, density and temperature
   */
  virtual void computeDimensionalPressDensTemp(CFreal& pDim, CFreal& rhoDim, CFreal& TDim, CFuint flxIdx);
  
protected:
  
  /// storage of the distance to the wall
  Framework::DataSocketSink<CFreal> socket_wallDistance;

  /// diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffVar;
  
  /// handle to qrad
  Framework::DataHandle<CFreal> m_qradFluxWall;  
    
  // flag telling if there is radiation coupling
  bool m_hasRadiationCoupling;

  /// arrray of gradients
  std::vector<RealVector*> m_gradients;
  
  //temporary value of density
  CFreal m_rhoWall;

  //temporary value of dyn viscosity
  CFreal m_muWall;

  //y+ value
  CFreal m_yPlus;
  
  //tau at the wall in 2D
  CFreal m_tau;
  
  // radiation heat flux
  CFreal m_heatFluxRad;
  
  //tau at the wall in 3D
  RealMatrix m_tau3D;
  
  //skin friction tensor in 3D
  RealMatrix m_Cf3D;
  
  // ID to identify the stanton number formula to use
  CFuint m_stantonNumID;
  
}; // end of class NavierStokesSkinFrictionHeatFluxFR

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxFR_hh
