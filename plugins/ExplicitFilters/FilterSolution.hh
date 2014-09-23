#ifndef COOLFluiD_Numerics_ExplicitFilters_FilterSolution_hh
#define COOLFluiD_Numerics_ExplicitFilters_FilterSolution_hh

//////////////////////////////////////////////////////////////////////////////

#include "FilterData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/SubSystemStatus.hh"
#include "MathTools/MathConsts.hh"
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
class FilterSolution : public FilterCom {
public:

  /**
   * Constructor.
   */
  explicit FilterSolution(const std::string& name) : 
		FilterCom(name),
  	socket_states("states"),
  	m_processRate()
  {
    addConfigOptionsTo(this);    // by default the data processing is run once -> processRate=infinity
    m_processRate = std::numeric_limits<CFuint>::max();
    setParameter("ProcessRate",&m_processRate);
  }

  /**
   * Destructor.
   */
  ~FilterSolution()
  {
  }
  
  static void defineConfigOptions(Config::OptionList& options)
   {
     options.addConfigOption< CFuint >("ProcessRate","Rate to process the data.");

   }

   void configure ( Config::ConfigArgs& args )
   {
     FilterCom::configure(args);
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

  CFuint m_processRate;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_FilterSolution_hh

