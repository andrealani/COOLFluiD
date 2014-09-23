#ifndef COOLFluiD_Physics_NEQ_NavierStokesNEQVarSet_hh
#define COOLFluiD_Physics_NEQ_NavierStokesNEQVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model for primitive
   * variables and chemical NEQ
   *
   * @author Andrea Lani
   * @author Janos Molnar
   */
template <typename BASE>
class NavierStokesNEQVarSet : public BASE {
public: // classes

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerNEQTerm;

  /**
   * Constructor
   * @see NavierStokes
   */
  NavierStokesNEQVarSet(const std::string& name,
			Common::SafePtr<Framework::PhysicalModelImpl> model);


  /**
   * Default destructor
   */

  virtual ~NavierStokesNEQVarSet();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  virtual void setComposition(const RealVector& state,
			      const bool isPerturb,
			      const CFuint iVar) = 0;

  /**
   * Set the quantities needed to compute gradients (pressure,
   * velocity, etc.) starting from the states
   */
  virtual void setGradientVars(const std::vector<RealVector*>& states,
			       RealMatrix& values,
			       const CFuint stateSize) = 0;

  /**
   * Get the adimensional dynamic viscosity
   */
  virtual CFreal getDynViscosity(const RealVector& state,
				 const std::vector<RealVector*>& gradients) = 0;

  /**
   * Get the adimensional density
   */
  virtual CFreal getDensity(const RealVector& state) = 0;

  /**
   * Get the adimensional thermal conductivity
   */
  virtual CFreal getThermConductivity(const RealVector& state,
				      const CFreal& dynViscosity) 
  {
    throw Common::NotImplementedException(FromHere(), "NavierStokesNEQVarSet::getThermConductivity()");
  }

  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius) = 0;

  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius) = 0;
  
protected:
  
  
  
  
  /// freeze enthalpies 
  void freezeEnthalpies(bool flag)
  {
    _freezeEnthalpies = flag; 
  }
  
protected:

  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// convective model
  Common::SafePtr<EulerNEQTerm> _eulerModel;

  /// freeze enthalpies flag
  bool _freezeEnthalpies;
  
  /// array for the mass composition
  RealVector _tempY;

  /// array to store the mass fractions of species
  RealVector _ys;

  /// matrix to store the diffusion velocities of species multiplied by the
  /// species densities
  RealVector _rhoUdiff;

  /// array to store the total enthalpies per unit mass of species
  RealVector _hsTot;

  /// molecules indices
  std::vector<CFuint> _moleculeIdx;

  /// array to store the vibrational enthalpies per unit mass of species
  RealVector _hsVib;

  /// array to store the electronic enthalpies per unit mass of species
  RealVector _hsEl;

  /// array to store the normal concentration gradients of species
  RealVector _normConcGradients;

  /// matrix to store the diffusion velocities of species multiplied by the
  /// species densities (backup)
  RealVector _rhoUdiffBkp;

  /// array to store the total enthalpies per unit mass of species (backup)
  RealVector _hsTotBkp;

  /// array to store the vibrational enthalpies per unit mass of species
  RealVector _hsVibBkp;

  /// array to store the electronic enthalpies per unit mass of species
  RealVector _hsElBkp;

  /// array to store the normal concentration gradients of species (backup)
  RealVector _normConcGradientsBkp;

  /// variable IDs (in the state) of the mass fractions
  std::vector<CFuint> _yID;

  /// variable IDs (in the state) of the vibrational temperatures
  std::vector<CFuint> _TvibID;

  /// average value of pressure
  CFreal _avP;
    
  /// vibrational temperatures
  RealVector _avTvib;

  /// dimensional vibrational temperature
  RealVector _avTvibdim;

  /// vibrational lambda
  RealVector _lambdaVib;

  /// Translational Rotational lambda
  CFreal _lambdaTR;

  /// non-dimensional thermal conductivity coefficient
  RealVector _thermCondCoeffVib;

}; // end of class NavierStokesNEQVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_NavierStokesNEQVarSet_hh
