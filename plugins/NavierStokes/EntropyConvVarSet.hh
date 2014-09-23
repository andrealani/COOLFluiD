#ifndef COOLFluiD_Physics_NavierStokes_EntropyConvVarSet_hh
#define COOLFluiD_Physics_NavierStokes_EntropyConvVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class implements an additioonal entropy equation to be plugged in 
  * an existing physical model
  *
  * @author Andrea Lani
  */
template <class BASE>
class EntropyConvVarSet : public BASE {
  
public: // classes
  
  /**
   * Constructor
   */
  EntropyConvVarSet(Common::SafePtr<Framework::BaseTerm> term);
  
  /**
   * Default destructor
   */
  virtual ~EntropyConvVarSet();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();
  
  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
				    RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                             RealVector& result);
  
  /**
   * Set other adimensional values for useful physical quantities
   */
  virtual void setDimensionalValuesPlusExtraValues
  (const Framework::State& state, RealVector& result,
   RealVector& extra);
  
protected:
  
  /**
   * Set enthalpy, energy, sound speed taking into account
   * thermal non-equilibrium effects
   */
  virtual void setEnthalpyEnergySoundSpeed(const Framework::State& state,
                                           RealVector& data,
                                           CFdouble& tempDim,
                                           CFdouble& pdim,
                                           CFdouble& rhodim);
    
}; // end of class EntropyConvVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "EntropyConvVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_EntropyConvVarSet_hh
