#ifndef COOLFluiD_Numerics_LUSGSMethod_LUFactorizationPivot_hh
#define COOLFluiD_Numerics_LUSGSMethod_LUFactorizationPivot_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "LUSGSMethod/LUFactorization.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class computes the LU factorization of the block diagonal Jacobian matrices, with pivoting.
   * @author Matteo Parsani
   * @author Kris Van den Abeele
   */
class LUFactorizationPivot : public LUFactorization {
public:

  /**
   * Constructor.
   */
  explicit LUFactorizationPivot(std::string name);

  /**
   * Destructor.
   */
  ~LUFactorizationPivot() {}

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks.
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: // functions

  void swapRows(const CFuint row1, const CFuint row2, RealMatrix& matrix);

protected: //Data

  /// socket for the pivot element of the LU factorization
  Framework::DataSocketSink< std::vector< CFuint > > socket_pivotLUFactorization;

}; // class LUFactorizationPivot

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_LUFactorizationPivot_hh

