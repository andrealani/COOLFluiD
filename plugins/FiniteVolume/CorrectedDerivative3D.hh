#ifndef COOLFluiD_Numerics_FiniteVolume_CorrectedDerivative3D_hh
#define COOLFluiD_Numerics_FiniteVolume_CorrectedDerivative3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "CorrectedDerivative2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes derivatives using a diamond volume in 3D
 *
 * @author Andrea Lani
 *
 */
class CorrectedDerivative3D : public CorrectedDerivative2D {
public:

  /**
   * Constructor
   */
  CorrectedDerivative3D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CorrectedDerivative3D();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:

  /*
   * Compute the inner gradients
   */
  virtual void computeInnerGradients(const RealMatrix& values,
				     std::vector<RealVector*>& gradients);

  /*
   * Compute the boundary gradients
   */
  virtual void computeBoundaryGradients(const RealMatrix& values,
					std::vector<RealVector*>& gradients);

  /*
   * Compute 1 - n x n
   */
  virtual void computeOEminusNN(const RealVector& n);
  
private: // data

  /// storage for uZ
  Framework::DataSocketSink<CFreal> socket_uZ;
  
  /// temporary array for RHS in LS algorithm
  RealVector _lf3;

}; // end of class CorrectedDerivative3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CorrectedDerivative3D_hh
