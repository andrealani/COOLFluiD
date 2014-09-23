#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCat_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCat_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSrvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
    class CatalycityModel;
  }
  
  namespace MathTools {
    class MatrixInverter;
  }

  namespace Numerics {
    
    namespace FiniteVolume {
            
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the Isothermal wall
   * BC in case of update variables including velocity components and temperature
   * [p v T y_i] or [rho_i v T] etc. Local Thermodynamic equilibrium is assumed 
   * at the wall, meaning that the species composition and the density at the wall
   * must be the equilibrium ones
   *  
   * @author Andrea Lani
   *
   */
template <class MODEL>
class NoSlipWallIsothermalNSrvtCat : public NoSlipWallIsothermalNSrvt<MODEL> {

public: 
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  NoSlipWallIsothermalNSrvtCat(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~NoSlipWallIsothermalNSrvtCat();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
protected:
  
  /**
   * Apply boundary condition on the given ghost state
   * @param innerState state vector associated to the internal cell
   * @param ghostState ghost state vector associated to the ghost cell
   * @param 
   */
  virtual void setGhostStateImpl(const Framework::State& innerState, 
				 Framework::State& ghostState);

private:
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// pointer to the convective term
  Common::SafePtr<MODEL> m_model;
 
  /// Numerical jacobian calculator
  std::auto_ptr<Framework::NumericalJacobian> m_numJacob;
  
  /// temporary data for holding the matrix inverter
  std::auto_ptr<MathTools::MatrixInverter> m_inverter;
  
  /// electron-vibrational temperatures
  RealVector m_tVec;
  
  /// mass species fractions in the internal cell
  RealVector m_ysIn;
  
  /// mass species fractions in the ghost cell
  RealVector m_ysG;
  
  /// mass species fractions at then wall
  RealVector m_ysWall;
  
  /// mass diffusion fluxes on the wall
  RealVector m_rhouDiff;
  
  /// mass production by catalycity on the wall
  RealVector m_catProduction;
  
  /// mass fractions gradients on the wall
  RealVector m_ysGradients;
  
  /// array of species mass fractions
  RealVector m_ysSol;

  /// residual 
  RealVector m_RHS;
  
  /// perturbed residual 
  RealVector m_pertRHS;
  
  /// residual difference 
  RealVector m_diffRHS;
  
  /// residual jacobian 
  RealMatrix m_dFdY;
  
  /// inverse of residual jacobian 
  RealMatrix m_invdFdY;
  
  /// catalycity model
  Common::SelfRegistPtr<Framework::CatalycityModel> m_catModel;
  
  /// name of the catalycity model
  std::string m_catModelName;
  
}; // end of class NoSlipWallIsothermalNSPvt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NoSlipWallIsothermalNSrvtCat.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallIsothermalNSrvtCat_hh
