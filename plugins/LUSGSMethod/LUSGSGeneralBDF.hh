#ifndef COOLFluiD_Numerics_LUSGSMethod_LUSGSGeneralBDF_hh
#define COOLFluiD_Numerics_LUSGSMethod_LUSGSGeneralBDF_hh

//////////////////////////////////////////////////////////////////////////////

#include "LUSGSMethod/LUSGSIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a ConvergenceMethod that implements the (nonlinear)
 * LU-SGS method applied to a general BDF time marching scheme. In this case the
 * coefficients of the time marching scheme are computed by the space method.
 *
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class LUSGSGeneralBDF : public LUSGSIterator {
public:

  /**
   * Default constructor without arguments.
   *
   * @param name missing documentation
   */
  explicit LUSGSGeneralBDF(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LUSGSGeneralBDF();

protected: // abstract interface implementations

  /**
   * Take one timestep.
   * @see ConvergenceMethod::takeStep()
   */
  virtual void takeStepImpl();

  /// Perform the prepare phase before any iteration
  virtual void prepare ();

}; // class LUSGSGeneralBDF

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_hh
