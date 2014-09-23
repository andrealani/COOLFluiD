#ifndef COOLFluiD_Physics_NavierStokes_Euler2DConsToCharInCons_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DConsToCharInCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

      class EulerTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to Roe to consistent conservative variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DConsToCharInCons : public Framework::VarSetMatrixTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DConsToCharInCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DConsToCharInCons();

  /**
   * Set the transformation matrix from a given state
   */
  void setMatrix(const RealVector& state);

private:

  /**
   * Set the flag telling if the transformation is an identity one
   * @pre this method must be called during set up
   */
  bool getIsIdentityTransformation() const
  {
    return false;
  }

private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

  /// temporary matrix
  RealMatrix _temp;

  /// matrix inverter size 4
  MathTools::MatrixInverterT<4> m_inverter4;

}; // end of class Euler2DConsToCharInCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DConsToCharInCons_hh
