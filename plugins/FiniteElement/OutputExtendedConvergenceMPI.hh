#ifndef COOLFluiD_Numerics_FiniteElement_OutputExtendedConvergenceMPI_hh
#define COOLFluiD_Numerics_FiniteElement_OutputExtendedConvergenceMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class outputs additional values at each iteration
 *
 * @author Thomas Wuilbaut
 */

class OutputExtendedConvergenceMPI : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  OutputExtendedConvergenceMPI(const std::string& name);

  /**
   * Default destructor
   */
  ~OutputExtendedConvergenceMPI();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

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
   * Execute on all dofs
   */
  void execute();

  /**
   * Compute the aerodynamic coefficients and write them to file
   */
  void computeExtraValues();

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFile();

  /**
   * Write the extra values to file
   */
  void updateOutputFile();

private:

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// Storage for choosing when to save the aero coef file
  CFuint _saveRate;

  /// Output File
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_outFile;

  /// Name of Output File
  std::string _nameOutputFile;

  RealVector _maxValues;

  RealVector _minValues;

  RealVector _averageState;

  RealVector _maxValuesGlobal;

  RealVector _minValuesGlobal;

  RealVector _averageStateGlobal;

}; // end of class OutputExtendedConvergenceMPI

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_OutputExtendedConvergenceMPI_hh
