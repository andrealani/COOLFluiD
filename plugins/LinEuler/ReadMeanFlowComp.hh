#ifndef COOLFluiD_Physics_LinEuler_ReadMeanFlowComp_hh
#define COOLFluiD_Physics_LinEuler_ReadMeanFlowComp_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/FunctionParser.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/VectorialFunction.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

   namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

///
/// This class create a socket for the mean flow variables and read the values from
/// a tecplot file (created by fluent for example)
/// @author Lilla Koloszar


class ReadMeanFlowComp : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  ReadMeanFlowComp(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ReadMeanFlowComp();

  /**
   * Configure the command
   */

  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();



  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();
   
  
  void writeMeanflowFile();
  
  void writeConnectivityFile();

protected: // functions

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

private: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket stores the data of the mean flow
  Framework::DataSocketSource<RealVector> socket_meanflow;

  /// the socket stores the states of the mean flow from Fluent
  Framework::DataSocketSource<RealVector> socket_fluent_states;

  /// the socket stores the coordinates of the mean flow from Fluent
  Framework::DataSocketSource<RealVector> socket_fluent_coords;

  /// string to hold the file name
  std::string m_file_name;

  /// number of nodes per element in Fluent
  std::valarray<CFuint> nbNodesPerElemFluent;
 
  
  bool _Interpolate;
  bool _Write;
  bool _WriteCT;
  
  ///connectivity table
  RealVector CTable;

}; // end of class ReadMeanFlowComp

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinEuler_ReadMeanFlowComp_hh
