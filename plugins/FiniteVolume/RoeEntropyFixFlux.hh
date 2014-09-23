#ifndef COOLFluiD_Numerics_FiniteVolume_RoeEntropyFixFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeEntropyFixFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "RoeFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux corresponding to the Euler
 * physical model 2D (in conservative variables)
 *
 * @author Andrea Lani
 *
 */
class RoeEntropyFixFlux : public RoeFlux {
public:

  /**
   * Constructor
   */
  RoeEntropyFixFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeEntropyFixFlux();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
private: // helper function
 
  /**
   * Compute the entropy-corrected lambda
   */
  void computeLambdaCorr(CFreal& lambdaCorr) const
  {
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i < nbEqs; ++i) {
      lambdaCorr = std::max(lambdaCorr,
       std::abs(_rightEvalues[i] - _leftEvalues[i]));
    }
  }
  
  /**
   * Set the abs of the eigen values
   */
  virtual void setAbsEigenValues();

private:
  
  /// ID of the entropy correction type
  CFuint _entropyFixID;
  
}; // end of class RoeEntropyFixFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeEntropyFixFlux_hh
