#ifndef COOLFluiD_Numerics_LUSGSMethod_LUSGSBDF2_hh
#define COOLFluiD_Numerics_LUSGSMethod_LUSGSBDF2_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>
using std::ofstream;

#include "LUSGSMethod/LUSGSIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a ConvergenceMethod that implements the (nonlinear)
 * LU-SGS method applied to a BDF2 time marching scheme.
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class LUSGSBDF2 : public LUSGSIterator {
public:

  /**
   * Default constructor without arguments.
   *
   * @param name missing documentation
   */
  explicit LUSGSBDF2(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LUSGSBDF2();

protected: // abstract interface implementations

  /**
   * Take one timestep.
   * @see ConvergenceMethod::takeStep()
   */
  virtual void takeStepImpl();

  /// Perform the prepare phase before any iteration
  virtual void prepare ();

protected: //data

  /// number of equations
  CFuint m_nbrEqs;

  /// temporay index iterator
  CFuint m_var_itr;

public: //data
  /// Name of Solution File where to write
  boost::filesystem::path fpath;

  /// Name of Convergence File where to write the convergence history.
  std::string m_nameConvergenceFile;

}; // class LUSGSBDF2

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_hh
