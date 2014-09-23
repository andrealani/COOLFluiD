#ifndef COOLFluiD_Numerics_LESDataProcessing_GradientComputer_hh
#define COOLFluiD_Numerics_LESDataProcessing_GradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalModel.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "LESProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/BaseDataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {
		  
		  class LESProcessingData;

//////////////////////////////////////////////////////////////////////////////

/** 
 * This class offers a basic interface for different gradient computer strategies.
 * 
 * @author Willem Deconinck
 */
class GradientComputer : public Framework::MethodStrategy<LESProcessingData> {

public: // functions

	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::BaseMethodStrategyProvider
  <LESProcessingData,GradientComputer > PROVIDER;

  /** 
   * Constructor
   */
  GradientComputer(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~GradientComputer();

  /** 
   * Set private data that will be used during the computation
   */
  virtual void setup();

  /** 
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "GradientComputer";
  }

  virtual void compute(std::vector<RealVector>& gradients, const CFuint& iCell) = 0 ;

  virtual CFreal getVolume(const CFuint& iCell) = 0 ;
  
  virtual CFreal getVolumeAdim(const CFuint& iCell) = 0;

//////////////////////////////////////////////////////////////////////////////


 
 	/**
	 * needsSockets()
	 * @return	a vector of SafePtr with the DataSockets needed as sinks
	 */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets()
    {
      CFAUTOTRACE;
      std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > result = m_sockets.getAllSinkSockets();
      std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > globalSockets = m_globalSockets.getAllSinkSockets();
      for(CFuint i=0; i<globalSockets.size(); ++i) {
        result.push_back(globalSockets[i]);
      }
      return result;
    }

  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets()
    {
      std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > result = m_sockets.getAllSourceSockets();
      std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > globalSockets = m_globalSockets.getAllSourceSockets();
      for(CFuint i=0; i<globalSockets.size(); ++i) {
        result.push_back(globalSockets[i]);
      }
      return result;
    }


  protected:

    /// Dynamic data sockets
    Framework::DynamicDataSocketSet<>                  m_sockets;
    Framework::DynamicDataSocketSet<Framework::GLOBAL> m_globalSockets;



}; // end of class GradientComputer

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_GradientComputer_hh
