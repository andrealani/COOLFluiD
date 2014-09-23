#ifndef COOLFluiD_Numerics_AnalyticalEE_ComputeIntegralError_hh
#define COOLFluiD_Numerics_AnalyticalEE_ComputeIntegralError_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/Node.hh"
#include "Framework/StdTrsGeoBuilder.hh"

#include "AnalyticalEE/AnalyticEEData.hh"
#include "Framework/ConvectiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a MethodCommand that computes the solution error
/// based on a user-defined analytic function.
/// @author Tiago Quintino
/// @author Martin Vymazal
class ComputeIntegralError : public AnalyticEECom {
public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit ComputeIntegralError(std::string name);

  /// Destructor.
  ~ComputeIntegralError();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Unset up private data and data of the aggregated classes
  /// in this command
  void unsetup();

  /// Configure the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

protected: // methods

  /// Compute error at quad point in reference space
  /// This method makes class ComputeIntegralError abstract base class

  virtual void integralErrorInCell(Framework::GeometricEntity *const cell, RealVector& error) = 0;

  /// Open the output file, write the header and close the file handle
  void prepareOutputFile();

  /// Open the output file
  /// @param bool indicates if reopen should be done with append
  void openFile(bool append);

  /// Compute the error at quadrature point


protected: // data

  ///data for numerical quadrature
  static CFreal qdWeightsTriag[7];
  static CFreal xQdPtsTriag[7];
  static CFreal yQdPtsTriag[7];
  static CFuint nbQdPts;
  static CFreal a1,b1,a2,b2,w1,w2,w3;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// the socket to the data handle of the nStatesProxy
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL > socket_states;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the variables
  std::vector<std::string> m_vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// Number of spatial dimensions
  CFuint m_dim;

  /// Number of equations
  CFuint m_nbEqs;

  /// vector to store coordinates of integration point in reference space
  RealVector m_mappedCoord;

  /// vector to store coordinates of integration point in physical space
  RealVector m_physicalCoord;

  /// vector to store the temporary value of the numerical solution
  RealVector m_numerical;

  /// vector to store the temporary value of the analytical solution
  RealVector m_exact;

  /// vector to store the temporary value of the error
  RealVector m_error;

  /// vector to store the temporary value of the L2_norm
  RealVector m_L2;

/// vector to store the temporary value of the logarithm of the L2_norm
  RealVector m_result;

/// temporary state to store numerical solution at quadrature point inside the element
  Framework::State* m_qdState;

  /// Output file for recording the evolution of the Error norm
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_file;

  /// Name of output file
  std::string m_nameOutputFile;

  /// Builder for standard TRS GeometricEntities
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder >  m_stdTrsGeoBuilder;

}; // class ComputeIntegralError

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AnalyticalEE_ComputeDiscreteError_hh

