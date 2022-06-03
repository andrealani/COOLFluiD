#ifndef COOLFluiD_FluxReconstructionMethod_ComputeFieldFromPotential_hh
#define COOLFluiD_FluxReconstructionMethod_ComputeFieldFromPotential_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a preprocessing command that computes the field from a potential in another subsystem
 *
 * @author Ray Vandenhoeck
 *
 */
class ComputeFieldFromPotential : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit ComputeFieldFromPotential(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ComputeFieldFromPotential();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions
  
  /**
   * Compute the flx pnt states
   */
  void computeFlxPntStates(std::vector< RealVector > states, std::vector< RealVector >& statesFlxPnt);

protected: // data

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// reconstruction coefficients for the flux points
  Common::SafePtr< RealMatrix > m_flxPntsRecCoefs;

  /// indexes of all the flux points
  Common::SafePtr< std::vector< CFuint > > m_allFlxPntIdxs;

  /// mapped coordinates of the solution points
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// dimensionality of the physical model
  CFuint m_dim;

  /// number of flux points for current element type
  CFuint m_nbrFlxPnts;

  /// number of solution points for current element type
  CFuint m_nbrSolPnts;
  
  /// maximum number of flux points in a cell
  CFuint m_maxNbrFlxPnts;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtNodes;
  
  /// coefs to compute the derivative of the states in the nodes
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtNodes;
  
  /// order of FR method
  CFuint m_order;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  /// extrapolated states in the flux points of the cell
  std::vector< RealVector > m_cellStatesFlxPnt;
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_otherUX;
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_otherUY;

  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_otherUZ;

  /// socket for nodes
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_otherStates;
  
  /// socket for state's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// storage of the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

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
  
  /// flag to use PFSS B as initialization
  bool m_usePFSSBInit;

}; // class ComputeFieldFromPotential

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ComputeFieldFromPotential_hh
