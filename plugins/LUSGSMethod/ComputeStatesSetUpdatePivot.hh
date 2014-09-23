#ifndef COOLFluiD_Numerics_LUSGSMethod_ComputeStatesSetUpdatePivot_hh
#define COOLFluiD_Numerics_LUSGSMethod_ComputeStatesSetUpdatePivot_hh

//////////////////////////////////////////////////////////////////////////////

#include "LUSGSMethod/ComputeStatesSetUpdate.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class updates the the states set solving two traingular after the LU factorization with pivoting.
   * @author Matteo Parsani
   * @author Kris Van den Abeele
   */
class ComputeStatesSetUpdatePivot : public ComputeStatesSetUpdate {
public:

  /**
   * Constructor.
   */
  explicit ComputeStatesSetUpdatePivot(std::string name);

  /**
   * Destructor.
   */
  ~ComputeStatesSetUpdatePivot() {}

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

protected:

  /// socket for the pivot element of the LU factorization
  Framework::DataSocketSink< std::vector< CFuint > > socket_pivotLUFactorization;

  /// pointer to the auxiliary rhs variable
  Common::SafePtr< RealVector > m_resAux;

}; // class ComputeStatesSetUpdate

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_ComputeStatesSetUpdatePivot_hh

