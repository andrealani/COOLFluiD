#ifndef COOLFluiD_FluxReconstructionMethod_BasePhysicality_hh
#define COOLFluiD_FluxReconstructionMethod_BasePhysicality_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command to check and enforce the physicality of the solution
 *
 * @author Ray Vandenhoeck
 *
 */
class BasePhysicality : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit BasePhysicality(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~BasePhysicality();

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

protected: // functions
  
  /**
   * Enforce physicality of the state
   */
  virtual void enforcePhysicality() = 0;
  
  /**
   * Check if the states are physical
   */
  virtual bool checkPhysicality() = 0;
  
  /**
   * Compute the flx pnt states
   */
  void computeFlxPntStates(std::vector< RealVector >& statesFlxPnt);

protected: // data

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// mapped coordinates of the solution points
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// dimensionality of the physical model
  CFuint m_dim;
  
  /// number of solution points for current element type
  CFuint m_nbrSolPnts;
  
  /// maximum number of flux points in a cell
  CFuint m_maxNbrFlxPnts;
  
  /// order of FR method
  CFuint m_order;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  /// extrapolated states in the flux points of the cell
  std::vector< RealVector > m_cellStatesFlxPnt;

}; // class BasePhysicality

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BasePhysicality_hh
