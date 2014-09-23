#ifndef COOLFluiD_Numerics_ExplicitFilters_StencilComputer_hh
#define COOLFluiD_Numerics_ExplicitFilters_StencilComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StencilComputerStrategy.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "ExplicitFilters/FilterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {
      
      class FilterData;
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class provides an abstract interface for functors
   * computing the stencil for a specific kind of polynomial
   * reconstruction in ExplicitFilters
   *
   * @author Willem Deconinck
   */
class StencilComputer : public Framework::StencilComputerStrategy<FilterData> {

public:

  typedef Framework::BaseMethodStrategyProvider<FilterData,StencilComputer> PROVIDER;
    
    // typedef Framework::MethodStrategyProvider<StencilComputer,
    //                                    FilterData, 
    //                                    StencilComputer,
    //                                    // Framework::StencilComputerStrategy< FilterData > , 
    //                                    ExplicitFiltersModule > PROVIDER
    //                                         // FilterData,
    //                                         //   Framework::StencilComputerStrategy<FilterData>,
    //                                         //   ExplicitFiltersModule> PROVIDER;

  /**
   * Returns the DataSockets that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
  providesSockets();
  
  /**
   * Returns the DataSockets that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /**
   * Defines the config options of this class
   * @param   options   config options of this class
   */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Configure the object
   * @param    args   configuration arguments
   */
  virtual void configure(Config::ConfigArgs& args );

  /**
   * Constructor
   */
   StencilComputer(const std::string& name);

  /**
   * Destructor
   */
  virtual ~StencilComputer();

  virtual void setup();


  virtual void prepareComputations();

  virtual void compute();
  virtual void computeInspected();
  virtual void computeStencil(const CFuint& iStencil);
  
  virtual void computeNeighbors(const CFuint& centreStateID, const CFuint& currStateID);

  virtual void postProcessStencil(const CFuint& centreStateID);

  void writeToFileStream(std::ofstream& fout);

  virtual void outputStencil();
  
  void computeWithLargerRadius(Common::SafePtr<std::vector<CFuint> > stencilsToCompute, CFreal enlargementFactor); 

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "StencilComputer";
  }
  
  bool mustPrecompute() { return m_precompute;}
  
  CFuint nbStencilsToCompute();
  CFuint nbInspectedStencilsToCompute();

protected: //data
  
  /// storage for the ghost states
  Framework::DataSocketSink < Framework::State*> socket_gstates;
  
private:
  
  bool m_precompute;
  bool m_prepared;
  bool m_computed;
  CFuint m_nbDistanceZones;
  CFuint m_nbBasicFilters;
  CFreal m_inclusionCriteria;
    
}; // end of class StencilComputer

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_StencilComputer_hh
