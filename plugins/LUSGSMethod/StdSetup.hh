#ifndef COOLFluiD_Numerics_LUSGSMethod_StdSetup_hh
#define COOLFluiD_Numerics_LUSGSMethod_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "LUSGSMethod/LUSGSIteratorData.hh"
#include "Framework/DataSocketSink.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to setup the MeshData.
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class StdSetup : public LUSGSIteratorCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdSetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSource<
      CFreal> socket_updateCoeff;

  /// socket for PastStates's
  Framework::DataSocketSource< Framework::State*> socket_pastStates;

  /// socket for current states set index
  Framework::DataSocketSource< CFint > socket_statesSetIdx;

  /// socket for bStatesNeighbors
  /// It is a list of the neighbor states for the boundary states.
  /// It will be useful to avoid very expensive jacobian matrix
  /// reallocations when applying strong boundary condition
  /// @note KVDA: don't think this is needed, at least not for SpectralFV/SpectralFD
//   Framework::DataSocketSource<
//                               std::valarray<Framework::State*> > socket_bStatesNeighbors;

  /// socket for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for the pivot element of the LU factorization
  Framework::DataSocketSource< std::vector< CFuint > > socket_pivotLUFactorization;

  /// dummy socket for rhs (not used)
  Framework::DataSocketSource< CFreal > socket_rhs;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_StdSetup_hh
