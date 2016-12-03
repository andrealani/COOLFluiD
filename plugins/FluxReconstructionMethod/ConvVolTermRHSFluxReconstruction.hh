#ifndef COOLFluiD_FluxReconstructionMethod_ConvVolTermRHSFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_ConvVolTermRHSFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/BaseVolTermComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic command that computes the volume terms for the
 * spectral finite difference schemes for convective terms
 *
 * @author Kris Van den Abeele
 *
 */
class ConvVolTermRHSFluxReconstruction : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit ConvVolTermRHSFluxReconstruction(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ConvVolTermRHSFluxReconstruction();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

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
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// set the data required to compute the volume term
  void setVolumeTermData();

  /// resize the variable m_resUpdates and m_gradUpdates corresponding to the current element type
  void resizeResAndGradUpdates();

  /// add the residual updates to the RHS
  void updateRHS();

  /// compute the volume term contribution to the gradients
  void computeGradients();

  /// add updates to gradients
  void addGradVolTermsAndDivideByJacobDet();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// Strategy that computes the volume terms in a cell
  Common::SafePtr< BaseVolTermComputer > m_volTermComputer;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// index of element type
  CFuint m_iElemType;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// updates to the residuals
  RealVector m_resUpdates;

  /// updates to the gradients
  std::vector< std::vector< RealVector > > m_gradUpdates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// number of dimensions in the physical model
  CFuint m_dim;

}; // class ConvVolTermRHSFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvVolTermRHSFluxReconstruction_hh
