#ifndef COOLFluiD_FluxReconstructionMethod_BCNoSlipWallrvt_hh
#define COOLFluiD_FluxReconstructionMethod_BCNoSlipWallrvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Physics {
    namespace NavierStokes {
      class EulerTerm;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a no-slip wall boundary condition for the 2D Euler/Navier-Stokes equations with rvt variables
 *
 * @author Ray Vandenhoeck
 * @author Firas Ben Ameur
 */
class BCNoSlipWallrvt : public BCStateComputer {

public:  // methods

  /// Constructor
  BCNoSlipWallrvt(const std::string& name);

  /// Destructor
  ~BCNoSlipWallrvt();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallrvt";
  }
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Set up private data and data
  void setup();

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                             std::vector< std::vector< RealVector* > >& ghostGrads,
                             const std::vector< RealVector >& normals,
                             const std::vector< RealVector >& coords);

protected: // data

  /// physical model (in conservative variables)
  Common::SafePtr< Framework::MultiScalarTerm< Physics::NavierStokes::EulerTerm > > m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;
  
  /// number of equations
  CFuint m_nbrEqs;
  
  /// physico-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// flag telling if the state has partial densities
  bool m_stateHasPartialDensities;
  
  /// number of species
  CFuint m_nbSpecies;
  
  /// number of vibrational temperatures
  CFuint m_nbTv;
  
  /// roto-translational and vibrational temperatures in the ghost state
  RealVector m_ghostTTvib;
  
  /// roto-translational and vibrational temperatures in the inner state
  RealVector m_innerTTvib;
  
  /// temperature var ID
  CFuint m_tempID;
  
  /// velocity IDs
  std::vector< CFuint > m_velocityIDs; 
  
  /// array of flags telling if a variable is a velocity component
  std::valarray<bool> m_isVelocityComp;
  
  /// wall static temperature
  CFreal m_wallT;

  /// iteration after which is changed to an isothermal wall BC
  CFuint m_changeToIsoT;

}; // class BCNoSlipWallrvt

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_BCNoSlipWallrvt_hh

