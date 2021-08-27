#ifndef COOLFluiD_Numerics_FiniteVolume_ComputeBField_hh
#define COOLFluiD_Numerics_FiniteVolume_ComputeBField_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
    namespace Numerics {
      
      namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class implements a pre-processing step during which the B field is read 
 * from file 
 *
 * @author Andrea Lani
 *
 */
class ComputeBField : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  ComputeBField(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ComputeBField();

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

  /**
   * Executes the mesh fitting
   */
  virtual void execute();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure(Config::ConfigArgs& args);

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:
  
  /**
   * Read the file defining the vectorial functions
   */
  virtual void readFunctionsFile();
  
private: //data
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// storage of face centroids
  Framework::DataSocketSink<CFreal> socket_faceCenters;
  
  /// socket for the Bfield in cells
  Framework::DataSocketSource<CFreal> socket_Bfield;

  /// socket for the Bfield in faces
  Framework::DataSocketSource<CFreal> socket_BfieldFaces;

  // dummy state vector
  RealVector m_input;
  
  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string> m_vars;
  
  /// name of the file where strings holding the functions are defined
  std::string m_functionsFileName;
  
  // the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;
  
}; // end of class ComputeBField

//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ComputeBField_hh
