#ifndef COOLFluiD_Numerics_FiniteElement_NewmarkComputeTimeResidual_hh
#define COOLFluiD_Numerics_FiniteElement_NewmarkComputeTimeResidual_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Common/CFMap.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to ComputeTimeResidual the MeshData.
 */
class NewmarkComputeTimeResidual : public FiniteElementMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit NewmarkComputeTimeResidual(const std::string& name);

  /**
   * Destructor.
   */
  ~NewmarkComputeTimeResidual();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Execute Processing actions
   */
  void executeOnTrs();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   *  Gets Alpha
   */
  CFreal getAlpha()
  {
  return _alpha;
  }

  /**
   *  Gets Gamma
   */
  CFreal getGamma()
  {
  return _gamma;
  }

private: // data

  /// map of LSSMatrix accumulators, one for each cell type
  Common::CFMap<CFuint,Framework::BlockAccumulator*> _mapAcc;

  /// socket for Rhs
  Framework::DataSocketSink<
                              CFreal> socket_rhs;

  /// socket for Rhs
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStates;

  /// socket for Rhs
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStatesD;

  /// socket for Rhs
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStatesD2;

  /// socket for flag for states on which strong BC is applied
  Framework::DataSocketSink<
                            std::vector<bool> > socket_appliedStrongBC;

  /// Storage of the Alpha Coefficient
  CFreal _alpha;

  /// Storage of the Gamma Coefficient
  CFreal _gamma;

  ///flag for using the row-sum mass matrix lumping
  bool _lumpMassMatrix;

  /// Temporary storage of the integration result
  RealMatrix _integResult;

}; // class NewmarkComputeTimeResidual

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NewmarkComputeTimeResidual_hh

