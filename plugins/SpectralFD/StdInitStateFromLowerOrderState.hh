#ifndef COOLFluiD_Numerics_SpectralFD_StdInitStateFromLowerOrderState_hh
#define COOLFluiD_Numerics_SpectralFD_StdInitStateFromLowerOrderState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

    class SpectralFDElementData;

//////////////////////////////////////////////////////////////////////////////

/**
 * An initalizing solution command for the spectral finite difference method,
 * which reads a lower order solution from a file and transfers it to the current order
 *
 * @author Kris Van den Abeele
 *
 */
class StdInitStateFromLowerOrderState : public SpectralFDMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit StdInitStateFromLowerOrderState(const std::string& name);

  /**
   * Destructor.
   */
  ~StdInitStateFromLowerOrderState();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * Execute Processing actions
   */
  void executeOnTrs();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Reads the CFmesh file header and performs some file checks
   */
  void readCFmeshFileHeaderData(std::ifstream& inputFile);

  /**
   * Computes the solution transformation matrix
   */
  void computeSolutionTransformationMatrix();

protected: // data

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> m_varSet;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// states in initialization points
  std::vector< Framework::State* > m_initPntsStates;

  /// input state
  Framework::State* m_inputState;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// dimensionality
  CFuint m_dim;

  /// coordinates of initialization point
  RealVector m_initPntCoords;

  /// name of the file containing the lower order solution
  std::string m_srcSolCFmeshFileName;

  /// polynomial order in the source CFmesh file
  CFuint m_srcSolPolyOrder;

  /// SpectralFDElementData object corresponding to solution polynomial order in input file
  SpectralFDElementData* m_srcSolSDData;

  /// transformation matrix from solution in input file to current solution
  RealMatrix m_solTransMatr;

  /// states in one cell in the input file
  std::vector< RealVector > m_cellSrcStates;
  
}; // class StdInitStateFromLowerOrderStateQuad

//////////////////////////////////////////////////////////////////////////////

    } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SpectralFD_StdInitStateFromLowerOrderState_hh

