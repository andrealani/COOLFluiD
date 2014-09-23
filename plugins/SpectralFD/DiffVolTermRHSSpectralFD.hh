#ifndef COOLFluiD_SpectralFD_DiffVolTermRHSSpectralFD_hh
#define COOLFluiD_SpectralFD_DiffVolTermRHSSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic command that computes the volume terms for the
 * spectral finite difference schemes for diffusive terms
 *
 * @author Kris Van den Abeele
 *
 */
class DiffVolTermRHSSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit DiffVolTermRHSSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DiffVolTermRHSSpectralFD();

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

  /// set the pointers to the gradients in the cell
  void setGradients();

  /// set the data required to compute the volume term
  void setVolumeTermData();

  /// resize the variables m_resUpdates and m_cellGrads corresponding to the current element type
  void resizeResUpdatesAndCellGrads();

  /// add the residual updates to the RHS
  void updateRHS();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// Strategy that computes the volume terms in a cell
  Common::SafePtr< BaseVolTermComputer > m_volTermComputer;

  /// index of element type
  CFuint m_iElemType;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// variable containing pointers to the gradients in a cell
  std::vector< std::vector< RealVector >* > m_cellGrads;

  /// updates to the residuals
  RealVector m_resUpdates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class DiffVolTermRHSSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_DiffVolTermRHSSpectralFD_hh
