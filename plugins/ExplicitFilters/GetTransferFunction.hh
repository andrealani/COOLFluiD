#ifndef COOLFluiD_Numerics_ExplicitFilters_GetTransferFunction_hh
#define COOLFluiD_Numerics_ExplicitFilters_GetTransferFunction_hh

//////////////////////////////////////////////////////////////////////////////

#include "FilterData.hh"
#include "Framework/SubSystemStatus.hh"

// #include "FilterCom.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

/**
 * A new Numerical Command 
 *
 * @author Willem Deconinck
 *
 */
class GetTransferFunction : public FilterCom {
public:

  /**
   * Constructor.
   */
  explicit GetTransferFunction(const std::string& name) 
  : FilterCom(name)
  {
    addConfigOptionsTo(this);
    // by default the data processing is run once -> processRate=infinity
    m_processRate = MathTools::MathConsts::CFrealMax();
    setParameter("ProcessRate",&m_processRate);
  }



  /**
   * Destructor.
   */
  ~GetTransferFunction()
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

  CFuint m_processRate;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_GetTransferFunction_hh

