#ifndef COOLFluiD_LES_LES2DVarSet_hh
#define COOLFluiD_LES_LES2DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace LES {

//////////////////////////////////////////////////////////////////////////////

/**
  * This class represents a 2D Large Eddy Simulation physical model
  *
  * @author Kris Van den Abeele
  * @author Ghader Ghorbaniasl
  */
class LES2DVarSet : public LESVarSet {
public: // classes

  /**
   * Constructor
   * @see LES2D
   */
  LES2DVarSet(const std::string& name,
              Common::SafePtr<Framework::PhysicalModelImpl> model) :
    LESVarSet(name, model)
  {
  }

  /**
   * Default destructor
   */
  virtual ~LES2DVarSet()
  {
  }

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);

  /**
   * Get the diffusive flux vector
   */
  virtual RealMatrix& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const CFreal& radius);

  /**
   * Get the heat flux
   */
  virtual CFreal getHeatFlux(const RealVector& state,
                             const std::vector<RealVector*>& gradients,
                             const RealVector& normal);

  /**
   * Get the axisymmetric source term
   */
  virtual void getAxiSourceTerm(const RealVector& physicalData,
                                const RealVector& state,
                                const std::vector<RealVector*>& gradients,
                                const CFreal& radius,
                                RealVector& source);

}; // end of class LES2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace LES

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_LES_LES2DVarSet_hh
