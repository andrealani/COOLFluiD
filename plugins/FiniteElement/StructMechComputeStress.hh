#ifndef COOLFluiD_Numerics_FiniteElement_StructMechComputeStress_hh
#define COOLFluiD_Numerics_FiniteElement_StructMechComputeStress_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Storage.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataProcessingData.hh"
#include "StructMech/StructMech2DDiffusiveDisp.hh"
#include "StructMech/StructMech2DDisp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Stresses
 * for StructMech2DDiffusiveDisp
 * These stresses are then outputted by the OutputFormatter
 *
 * @author Thomas Wuilbaut
 *
 */
class StructMechComputeStress : public Framework::DataProcessingCom {
public:

  /**
   * Constructor.
   */
  StructMechComputeStress(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMechComputeStress();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  void unsetup();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

private: //data

  /// Variable Set
  Common::SelfRegistPtr<Physics::StructMech::StructMech2DDiffusiveDisp> _varSet;

  /// storage for stress
  Framework::DataSocketSource<RealVector> socket_stress;

  /// storage for strain
  Framework::DataSocketSource<RealVector> socket_strain;

  /// storage for strain
  Framework::DataSocketSource<CFuint> socket_flagStates;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
}; // end of class StructMechComputeStress

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_StructMechComputeStress_hh
