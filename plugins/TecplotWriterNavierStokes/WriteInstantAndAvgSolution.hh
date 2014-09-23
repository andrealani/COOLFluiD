#ifndef COOLFluiD_IO_TecplotWriter_WriteInstantAndAvgSolution_hh
#define COOLFluiD_IO_TecplotWriter_WriteInstantAndAvgSolution_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/FileWriter.hh"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/VarSetTransformer.hh"

#include "TecplotWriter/WriteSolution.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action
 * to write the solution in Pvt variables to a Tecplot format file
 * for visualization. An additional file with the averaged
 * solution (for LES) until now is written.
 *
 * @note to be used only with Navier-Stokes
 * @warning Average solution file only implemented for ASCII.
 *
 * @author Kris Van den Abeele
 *
 */
class WriteInstantAndAvgSolution : public WriteSolution {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit WriteInstantAndAvgSolution(const std::string& name);

  /**
   * Destructor.
   */
  ~WriteInstantAndAvgSolution()
  {
  }

  /**
    * Set up private data
   */
  void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

protected:

  /**
   * Write the to the given file stream the MeshData.
   * @throw Common::FilesystemException
   */
  void writeToFileStream(std::ofstream& fout);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected:

  /// socket for averaged variables in nodes
  Framework::DataSocketSource<RealVector> socket_nodeAvgVars;

  /// Transformer from update to primitive Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_updateToPrimVar;

  /// Update variable set string
  /// @note it would be better if this could somehow be gotten from the spatial method instead of letting the user define it
  std::string m_updateVarStr;

  /// Pvt state variable names
  std::vector< std::string > m_primVarNames;

  /// extra averaged variable names
  std::vector< std::string > m_extraAvgVarNames;

  /// variable for nodal Pvt state
  Framework::State* m_nodalPrimState;

  /// variable for averaged values in a node
  RealVector m_nodeAvgVals;

  /// counter for the number of averaging steps
  CFuint m_avgStepsCounter;

  /// rate for writing the instantaneous and averaged solution to the file
  CFuint m_writeToFileRate;

  /// boolean telling whether the current instantaneous solution should be added to the average solution
  /// (this allows to omit the initial solution in the case of a restart)
  bool m_addCurrSolToAvgSol;

  /// name of the file for the averaged solution
  std::string m_avgSolFileNameStr;

  /// name of averaged solution file where to write
  boost::filesystem::path m_avgSolFileName;

}; // class WriteInstantAndAvgSolution

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_TecplotWriter_WriteInstantAndAvgSolution_hh

