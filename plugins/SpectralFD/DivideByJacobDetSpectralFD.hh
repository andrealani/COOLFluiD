#ifndef COOLFluiD_SpectralFD_DivideByJacobDetSpectralFD_hh
#define COOLFluiD_SpectralFD_DivideByJacobDetSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a SpectralFD command that multiplies the update coefficients with the
 * Jacobian determinant and divides by the cell volume (because they are not really the same...)
 * This command is only needed for steady computations (for unsteady, volumes are computed for each state,
 * and the volumes datasocket actually contains the Jacobian determinants)
 *
 * @author Kris Van den Abeele
 *
 */
class DivideByJacobDetSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit DivideByJacobDetSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DivideByJacobDetSpectralFD();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unsetup private data
   */
  void unsetup();

  /**
   * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected: // functions

  /**
   * Divides by jacobian determinant
   */
  void divideRHSByJacobDet();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class DivideByJacobDetSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_DivideByJacobDetSpectralFD_hh
