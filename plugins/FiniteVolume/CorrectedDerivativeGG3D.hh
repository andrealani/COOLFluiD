#ifndef COOLFluiD_Numerics_FiniteVolume_CorrectedDerivativeGG3D_hh
#define COOLFluiD_Numerics_FiniteVolume_CorrectedDerivativeGG3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/CorrectedDerivativeGG2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes derivatives using a diamond volume in 3D
 *
 * @author Andrea Lani
 *
 */
class CorrectedDerivativeGG3D : public CorrectedDerivativeGG2D {
public:

  /**
   * Constructor
   */
  CorrectedDerivativeGG3D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CorrectedDerivativeGG3D();

  /**
   * Set up the member data
   */
  virtual void setup();

protected:

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

  /// temporary array for RHS in LS algorithm
  RealVector _lf3;

}; // end of class CorrectedDerivativeGG3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CorrectedDerivativeGG3D_hh
