#ifndef COOLFluiD_Numerics_ExplicitFilters_Prepare_hh
#define COOLFluiD_Numerics_ExplicitFilters_Prepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "FilterData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/SubSystemStatus.hh"
// #include "FilterCom.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
	
	namespace Framework {
		class State;
	}

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

/**
 * Prepare Command 
 *
 * This Command Prepares the Filtering before computations start.
 * It computes all stencils and weights
 * 
 * @author Willem Deconinck
 *
 */
class Prepare : public FilterCom {
public:

  /**
   * Constructor.
   */
  explicit Prepare(const std::string& name) : 
		FilterCom(name),
		socket_states("states"),
		socket_nodes("nodes")
{
    addConfigOptionsTo(this);
    // by default the data processing is run once -> processRate=infinity
    m_processRate = std::numeric_limits<CFuint>::max();
    setParameter("ProcessRate",&m_processRate);
    m_outputBadFilters = true;
    setParameter("OutputBadFilters",&m_outputBadFilters);
    
  }

  /**
   * Destructor.
   */
  ~Prepare()
  {
  }
  
  static void defineConfigOptions(Config::OptionList& options)
   {
     options.addConfigOption< CFuint >("ProcessRate","Rate to process the data.");
     options.addConfigOption< bool >("OutputBadFilters","Output bad filters (default=true)");
   }

   void configure ( Config::ConfigArgs& args )
   {
     CFLog(INFO, "Configuring Prepare");
     Framework::MethodCommand<FilterData>::configure(args);
     CFLog(INFO, "... done \n");
   }
   
   /**
    * Returns the DataSockets that this command provides as sources
    * @return a vector of SafePtr with the DataSockets
    */
   virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
   needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();
  
  void outputBadFilters();

private:
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
    
  void writeBadFiltersToFileStream(std::ofstream& fout, std::ofstream& lout);

  CFuint m_processRate;
  
  bool m_outputBadFilters;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_Prepare_hh

