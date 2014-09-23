#ifndef COOLFluiD_LES_LESVarSet_hh
#define COOLFluiD_LES_LESVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace LES {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Large Eddy Simulation physical model
/// @author Kris Van den Abeele
/// @author Ghader Ghorbaniasl
class LESVarSet : public Physics::NavierStokes::NavierStokesVarSet {
public: // classes

  /// Constructor
  /// @see NavierStokes
  LESVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  virtual ~LESVarSet();

  /**  Set up private data */
  virtual void setup();

  /// Set the composition
  /// @pre this function has to be called before any other function
  ///     computing other physical quantities
  virtual void setComposition(const RealVector& state,
                              const bool isPerturb,
                              const CFuint iVar) {};

  /** Get the adimensional dynamic viscosity */
  virtual CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients) = 0;
  
  /** Get the dynamic SGS viscosity */
  virtual CFreal getDynSGSViscosity(const RealVector& state, 
                                    const std::vector<RealVector*>& gradients, 
                                    const CFreal& volume);
                                    
  /** Get the SGS kinetic energy */
  virtual CFreal getSGSKineticEnergy(const RealVector& state, 
                                    const std::vector<RealVector*>& gradients, 
                                    const CFreal& volume);                                  
                                    
  /** Get the adimensional density */
  virtual CFreal getDensity(const RealVector& state) = 0;

  /** Get the adimensional thermal conductivity */
  virtual CFreal getThermConductivity(const RealVector& state,
                                      const CFreal& dynViscosity) = 0;

  /// Puts the velocity on the passed variable
  /// @param velocity variable to store the velocity
  /// @pre setGradientState()
  void putVelocity(double* velocity);

  /// Puts the velocity gradients on the passed variable
  /// @param velGrads variable to store the velocity gradients
  void putVelGrads(double* velGrads);

  /// Puts the temperature on the passed variable
  /// @param temperature variable to store the temperature
  /// @pre setGradientState()
  void putTemperature(double* temperature);

  /// Puts the temperature gradient on the passed variable
  /// @param tempGrads variable to store the temperature gradient
  void putTempGrad(double* tempGrads);

  /// Puts the density on the passed variable
  /// @param density variable to store the density
  void putDensity(double* density);

  /// Puts the heat capacity at constant pressure on the passed variable
  /// @param cp variable to store the heat capacity
  void putCp(double* cp);

  /// Puts the volume on the passed variable
  /// @param vol variable to store the volume
  void putVolume( double* vol);
  
  /// Set the volume
  void setVolume(const CFreal& vol) { m_stateVolume = vol;}
  
protected:
  
  /// Get the volume
  CFreal getVolume() const { return m_stateVolume;}

protected: // variables

  /// dimensionality
  CFuint m_dimensionality;

  /// number of equations
  CFuint m_nbrEqs;

  /// number of equations minus one
  CFuint m_nbrEqsMin1;

  /// Euler model
  Common::SafePtr<Physics::NavierStokes::EulerTerm> m_eulerModel;

  /// current state
  Common::SafePtr< const RealVector > m_currState;

  /// current gradients
  Common::SafePtr< const std::vector< RealVector* > > m_currGrad;

private:
  /// volume associated to current state
  CFreal m_stateVolume;

}; // end of class LESVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace LES

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_LES_LESVarSet_hh
