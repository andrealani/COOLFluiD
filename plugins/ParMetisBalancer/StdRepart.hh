#ifndef COOLFluiD_IO_ParMetisBalancer_ReadCFmesh_hh
#define COOLFluiD_IO_ParMetisBalancer_ReadCFmesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "ParMetisBalancer/ParMetisBalancerData.hh"

#include "DataStorage.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace ParMetisBalancer {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class is a MethodCommand that dynamicaly balances the mesh
 *
 * @author
 *
 */
class StdRepart : public ParMetisBalancerCom {

public: // functions

  /**
   * Constructor
   */
  explicit StdRepart(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdRepart();

  /**
   * Configures the command.
   */
  void configure(Config::ConfigArgs& args);

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // helper functions

  /**
  * Setup info about the mesh
  * info about nodes, and whare thay belong
  */
  void setupDataStorage();
  
  /**
  * Builds global numbering for nodes, so it is understood by PARMetis
  */
  void globalNumNode();
  
  /**
  * Builds information about mesh conectivity
  * Stored in CSR format, such that ParMetis may understand
  * Works for simplexes
  */
  void setupCSR();

  /**
  * Parmetis ParMETIS_V3_AdaptiveRepart is called to calculate new mesch partitoning
  */
  void callParMetisAdaptiveRepart();

  /**
  * Processes comunicate their interface part flag
  */
  void UpdateInterfacePart2();

  /**
  * Make lists of elements for which certain action should be taken
  */
  void SelectCellsToSend();
  
  /**
  * Make lists of nodes for which certain action should be taken
  */
  void SelectNodesToSend();
  
  /**
  * Make MPI structs to be sent recived
  */
  void PrepereMPIcommStruct();
  
  /**
  * Do MPI comm
  */
  void MPICommunicate();
  
  /**
  * preper to update
  */
  void PrepereToUpdate();
  
  /**
  * preper to update - select nodes to delete
  */
  void SelectNodesToDelete();
  
//  /**
//  * Apply New Mesh - add/delete nodes/cells/states
//  */
//  void ApplyNewMesh();

  /**
  * dealocates memory after use
  */
  void DoClearMemory();

  /**
  *  My function to produce TecPlot output with data that interest me
  * TODO: this is to be removed later
  */
  template<CFuint DIM>
  void DoWriteTec(boost::filesystem::path& path, bool interface_view);
  
  /**
  *  My function to produce TecPlot output with data that interest me
  *  omits removed elements
  * TODO: this is to be removed later
  */
  template<CFuint DIM>
  void DoWriteTecNoRemoved(boost::filesystem::path& path);
  
  /**
  *  My function to produce TecPlot output with data that interest me
  *  prints data recived from others
  * TODO: this is to be removed later
  */
  vector<CFuint> toBeSavedNodesGlobalId; // container used to specify which cell of the orginal mesh are needed for new elements
  template<CFuint DIM>
  void DoWriteTecAfterSendRecive(boost::filesystem::path& path);

  /**
  *  Some other function is here...
  */

private: // data

  /// Repartition data storage
  DataStorage dataStorage;

  // ---------------------------------------

  /// the socket to the data handle of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  ///
  Common::SafePtr<Framework::TopologicalRegionSet> m_cells;

  ///
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > m_nodes;
  
  Framework::DataHandle < Framework::State*, Framework::GLOBAL > m_states;


}; // class ReadCFmesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParMetisBalancer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_ParMetisBalancer_ReadCFmesh_hh

