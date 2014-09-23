#ifndef COOLFluiD_Numerics_FiniteVolume_RoeFluxALE_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeFluxALE_hh

//////////////////////////////////////////////////////////////////////////////

#include "RoeFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux for Arbitrary Langragian Eulerian formulation
 * using the Cranck-Nicholson time integration and geometric parameter
 * averaging (see I.Lepot thesis, section 5.4)
 *
 * @author Thomas Wuilbaut 
 * @author Andrea Lani
 *
 */
class RoeFluxALE : public RoeFlux {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  RoeFluxALE(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeFluxALE();

  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Compute the flux in the current face
   */
  virtual void compute(RealVector& result);
  
protected:

  /**
   * Set the abs of the eigen values
   */
  virtual void setAbsEigenValues();
  
  /**
   * Compute the entropy-corrected lambda
   */
  void computeLambdaCorr(CFreal& lambdaCorr) const
  {
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i < nbEqs; ++i) {
      lambdaCorr = std::max(lambdaCorr, std::abs(_rightEvalues[i] - _leftEvalues[i]));
    }
  }
  
protected: //data
  
  /// storage of the past nodes
  Framework::DataSocketSink <Framework::Node*> socket_pastNodes;
  
  /// storage of the past states
  Framework::DataSocketSink <Framework::Node*> socket_futureNodes;
  
  /// Mesh Speed projected on the unit normal
  CFreal   _vgn;
  
  /// Mesh Speed
  RealVector   _meshSpeed;
  
  /// ID of the entropy correction type
  CFuint _entropyFixID;
  
}; // end of class RoeFluxALE

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeFluxALE_hh
