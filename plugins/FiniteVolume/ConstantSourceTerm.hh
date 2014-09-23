#ifndef COOLFluiD_Numerics_FiniteVolume_ConstantSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_ConstantSourceTerm_hh

//////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term of constant value
 *
 * @author Thomas Wuilbaut
 *
 */
class ConstantSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  ConstantSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~ConstantSourceTerm();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Compute the source term
   * @param element   this GeometricEntity may or may not be used 
   * @param cellID    ID can always be useful in case GeometricEntity is not created
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);
    
}; // end of class ConstantSourceTerm

//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ConstantSourceTerm_hh
