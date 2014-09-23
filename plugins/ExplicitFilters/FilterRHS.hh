#ifndef COOLFluiD_Numerics_ExplicitFilters_FilterRhs_hh
#define COOLFluiD_Numerics_ExplicitFilters_FilterRhs_hh

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
 * A new Numerical Command 
 *
 * @author Willem Deconinck
 *
 */
class FilterRhs : public FilterCom {
public:

  /**
   * Constructor.
   */
  explicit FilterRhs(const std::string& name) : 
		FilterCom(name),
  	socket_states("states"),
    socket_rhs("rhs")
  {
    addConfigOptionsTo(this);
    // by default the data processing is run once -> processRate=infinity
    m_processRate = std::numeric_limits<CFuint>::max();
    setParameter("ProcessRate",&m_processRate);
  }

  /**
   * Destructor.
   */
  ~FilterRhs()
  {
  }
  
  static void defineConfigOptions(Config::OptionList& options)
   {
     options.addConfigOption< CFuint >("ProcessRate","Rate to process the data.");
   }

   void configure ( Config::ConfigArgs& args )
   {
     Framework::MethodCommand<FilterData>::configure(args);
   }

  /**
   * Execute Processing actions
   */
  void execute();

	/**
	 * needsSockets()
	 * @return	a vector of SafePtr with the DataSockets needed as sinks
	 */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();


private:

  /// Handle to States
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
    
  /// handle to rhs
  Framework::DataSocketSink<CFreal> socket_rhs;


  CFuint m_processRate;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_FilterRhs_hh

