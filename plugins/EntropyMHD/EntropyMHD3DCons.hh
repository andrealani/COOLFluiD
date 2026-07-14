#ifndef COOLFluiD_Physics_EntropyMHD_EntropyMHD3DCons_hh
#define COOLFluiD_Physics_EntropyMHD_EntropyMHD3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "EntropyMHD3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * Conservative variable set for EntropyMHD 3D MHD with GLM projection.
 *
 * State: (rho, rhoU, rhoV, rhoW, Bx, By, Bz, sigma, phi)
 *
 * sigma = p*rho^(1-gamma) (entropy variable, same in Cons and Prim).
 *
 * @author Rayan Dhib
 */
class EntropyMHD3DCons : public EntropyMHD3DVarSet {
public:

  /**
   * Constructor
   */
  EntropyMHD3DCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~EntropyMHD3DCons();

  /**
   * Set up the private data
   */
  virtual void setup();

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const { return 9; }

  /**
   * Set the jacobians
   */
  void computeJacobians() {}

  /**
   * Compute the physical data from conservative state
   * @param state  conservative state vector
   * @param data   physical data vector (output)
   */
  void computePhysicalData(const Framework::State& state,
                           RealVector& data);

  /**
   * Compute the state from physical data
   * @param data   physical data vector
   * @param state  conservative state vector (output)
   */
  void computeStateFromPhysicalData(const RealVector& data,
                                    Framework::State& state);

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  void setDimensionalValues(const Framework::State& state,
                            RealVector& result)
  {
    // non-dimensional by default
    result = static_cast<const RealVector&>(state);
  }

  /**
   * Compute the perturbed states data
   */
  void computePerturbedPhysicalData(const Framework::State& state,
                                    const RealVector& bData,
                                    RealVector& data,
                                    CFuint iVar)
  {
    cf_assert(iVar < 9);
    computePhysicalData(state, data);
  }

  /// Set the IDs corresponding to the velocity components
  virtual void setStateVelocityIDs(std::vector<CFuint>& velIDs)
  {
    velIDs.resize(3);
    velIDs[XX] = 1;
    velIDs[YY] = 2;
    velIDs[ZZ] = 3;
  }

}; // end of class EntropyMHD3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_EntropyMHD_EntropyMHD3DCons_hh
