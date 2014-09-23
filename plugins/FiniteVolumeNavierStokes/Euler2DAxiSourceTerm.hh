#ifndef COOLFluiD_Numerics_FiniteVolume_Euler2DAxiSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_Euler2DAxiSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 2D for conservative
 * variables
 *
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class Euler2DAxiSourceTerm : public ComputeSourceTermFVMCC {

public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  Euler2DAxiSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~Euler2DAxiSourceTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
  }
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);
    
private: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;

  /// vector to store temporary result
  RealVector _temp;

  /// Euler physical data
  RealVector _physicalData;

  /// ID of the equation to set
  CFuint _sID;
  
}; // end of class Euler2DAxiSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Euler2DAxiSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Euler2DAxiSourceTerm_hh
