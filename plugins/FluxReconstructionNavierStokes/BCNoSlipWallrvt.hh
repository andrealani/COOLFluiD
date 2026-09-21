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
 * Ghost state, with u_i the interior state at the wall flux point (the trace of the cell polynomial):
 * velocity reversed, each temperature reflected about the wall temperature on an isothermal wall
 * (T_g = 2 T_wall - T_i, at least 0.01 T_wall) or copied before ChangeToIsoT, species densities copied
 * (partial pressures scaled by T_g / T_i when the state carries them), the electron entry and Te copied.
 * LegacyGhost = true restores the previous ghost (wall temperature, densities scaled by T_i / T_g).
 *
 * @author Ray Vandenhoeck
 * @author Firas Ben Ameur
 * @author Rayan Dhib
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

  /**
   * Sets the boundary values of the gradient variables: zero velocity, the wall
   * temperatures that the wall sets and, for the species and the adiabatic
   * temperatures, the gradient variables extrapolated to the flux points.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states: the wall temperatures that the wall sets, the
   * partial densities they give from the interior partial pressures, and the
   * average of the interior and ghost states for the velocity and the adiabatic
   * temperatures.
   */
  void computeBndStates(const std::vector< Framework::State* >& intStates,
                        const std::vector< Framework::State* >& ghostStates,
                        const std::vector< RealVector >& unitNormals,
                        const std::vector< RealVector >& flxPntCoords,
                        std::vector< RealVector* >& bndStates);

  /**
   * Sets the boundary gradients: the compact face gradients, without the normal
   * components of T and Tv before ChangeToIsoT and without the normal components
   * of the mass fractions when NonCatalytic is true.
   */
  void computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                       std::vector< std::vector< RealVector* > >& bndGrads,
                       const std::vector< RealVector* >& bndStates,
                       const std::vector< RealVector >& unitNormals,
                       const std::vector< RealVector >& flxPntCoords);

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

  /// how many non finite ghost states have already been reported
  CFuint m_nbBadGhostReported;
  
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

  /// use the previous ghost state instead of the reflected one
  bool m_legacyGhost;

  /// no species diffusion flux through the wall
  bool m_nonCatalytic;

}; // class BCNoSlipWallrvt

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_BCNoSlipWallrvt_hh

