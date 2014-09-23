#ifndef COOLFluiD_Numerics_FiniteVolumeMHD_iPic3DCoupler_hh
#define COOLFluiD_Numerics_FiniteVolumeMHD_iPic3DCoupler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Environment/FileHandlerOutput.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolumeMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the divergence of magnetic and electric fields
 *
 * @author Andrea Lani
 * @author Vyacheslav Olshevsky
 */
class iPic3DCoupler : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */

  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  iPic3DCoupler(const std::string& name);

  /**
   * Default destructor
   */
  ~iPic3DCoupler();

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

  /**
   * Execute on a set of dofs
   */
  void execute();

  /**
   * Configures this object with supplied arguments.
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
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
  providesSockets();
  
private: //data
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage for State's
  Framework::DataSocketSink < Framework::State* > socket_gstates;
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;
  
  /// socket for uZ values
  Framework::DataSocketSink<CFreal> socket_uZ;
    
  /// storage of curl B
  Framework::DataSocketSource < CFreal > socket_curlB;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;  
  
  /// update variable set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> m_updateVarSet;
  
  /// save rate
  CFuint  m_saveRate;
  
  /// Save each iteration to a different Tecplot file (with suffix _iter#).
  bool  m_appendIter;

  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  bool  m_appendTime;
  
  /// path to input file for iPic3D
  std::string m_pathToInputFile; 
  
}; // end of class iPic3DCoupler

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeMHD

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeMHD_iPic3DCoupler_hh
