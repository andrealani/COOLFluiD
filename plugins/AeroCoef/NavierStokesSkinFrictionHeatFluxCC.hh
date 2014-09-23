#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCC_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "AeroCoef/AeroForcesFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }

  namespace Numerics {
    namespace FiniteVolume {
      class DerivativeComputer;
    }
  }
  
  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the skin friction and the heat flux for NavierStokes
 * simulations with
 * @see CellCenterFVM
 *
 * @author Andrea Lani
 * @author Thomas Wuilbaut
 *
 */
class NavierStokesSkinFrictionHeatFluxCC : public AeroForcesFVMCC {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesSkinFrictionHeatFluxCC(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesSkinFrictionHeatFluxCC();
  
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
   * Compute dimensional pressure, density and temperature
   */
  virtual void computeDimensionalPressDensTemp(CFreal& pDim, CFreal& rhoDim, CFreal& TDim);
  
  /**
   * Compute the value of Y+
   */
  void computeYplus();
  
  /// Update the data to write
  virtual void updateWriteData();
  
  /**
   * Compute the value of tau at the wall
   */
  void computeTauWall();
  
protected:
  
  /// storage of the distance to the wall
  Framework::DataSocketSink<CFreal> socket_wallDistance;
  
  // derivative computer
  Common::SafePtr<Numerics::FiniteVolume::DerivativeComputer> _derivComputer;

  /// diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffVar;
  
  /// handle to qrad
  Framework::DataHandle<CFreal> _qradFluxWall;  
    
  // flag telling if there is radiation coupling
  bool _hasRadiationCoupling;
  
  // array of states involved in the flux computation
  std::vector<RealVector*> _states;
  
  // array of values (p, u, v, T, ...)
  RealMatrix _values;

  /// arrray of gradients
  std::vector<RealVector*> _gradients;
  
  //temporary value of density
  CFreal _rhoWall;

  //temporary value of dyn viscosity
  CFreal _muWall;

  //y+ value
  CFreal _yPlus;
  
  //tau at the wall in 2D
  CFreal _tau;
  
  // radiation heat flux
  CFreal _heatFluxRad;
  
  //tau at the wall in 3D
  RealMatrix _tau3D;
  
  //skin friction tensor in 3D
  RealMatrix m_Cf3D;
  
  // ID to identify the stanton number formula to use
  CFuint _stantonNumID;
  
}; // end of class NavierStokesSkinFrictionHeatFluxCCNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCC_hh
