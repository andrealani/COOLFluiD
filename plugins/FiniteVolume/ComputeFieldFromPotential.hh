#ifndef COOLFluiD_Numerics_FiniteVolume_ComputeFieldFromPotential_hh
#define COOLFluiD_Numerics_FiniteVolume_ComputeFieldFromPotential_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////

/**
 * This class computes a field by taking gradients of a potential field
 *
 * @author Andrea Lani
 *
 */
class ComputeFieldFromPotential : public CellCenterFVMCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  explicit ComputeFieldFromPotential(const std::string& name);
  
  /**
   * Default destructor
   */
  ~ComputeFieldFromPotential();

  /**
   * Configure
   */
  void configure ( Config::ConfigArgs& args );
  
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

  /**
   * Execute on a set of dofs
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  // std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  // storage of the past states
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_otherUX;
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_otherUY;

  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_otherUZ;

  /// socket for nodes
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_otherStates;

  /// flag telling whether to apply the post processing
  bool m_applyProcessing;
  
  /// IDs of the state variables to assign to the newly computed field
  std::vector<CFuint> m_variableIDs;

  /// name of the other namespace (providing the potential)
  std::string m_otherNamespace;

  /// radius corresponding to the internal boundary between donor and current grids
  CFreal m_interRadius;
  
  /// distance within which points in the smaller mesh are selected
  CFreal m_deltaSelection;
  
}; // end of class ComputeFieldFromPotential
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ComputeFieldFromPotential_hh



