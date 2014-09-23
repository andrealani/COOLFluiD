#ifndef COOLFluiD_Numerics_LESDataProcessing_TurbulenceFunction_hh
#define COOLFluiD_Numerics_LESDataProcessing_TurbulenceFunction_hh

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
 * This class offers a basic interface for different turbulence function strategies.
 * 
 * @author Willem Deconinck
 */
class TurbulenceFunction : public Framework::MethodStrategy<LESProcessingData> {

public: // functions
  
	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::BaseMethodStrategyProvider
  <LESProcessingData,TurbulenceFunction > PROVIDER;

  /** 
   * Constructor
   */
  TurbulenceFunction(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~TurbulenceFunction();

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
    return "TurbulenceFunction";
  }
  
  virtual void setSocketName() = 0;

  virtual void compute(const RealVector& states, std::vector<RealVector>& gradients, const CFuint& iCell) = 0 ;
  
  virtual void computeAverage(const CFreal& oldWeight, const CFreal& newWeight, const RealVector& states, std::vector<RealVector>& gradients, const CFuint& iCell);
  
  virtual bool extrapolate() 
  {
    return m_extrapolate;
  }
  
  virtual CFuint getNbStates()
  {
    return m_nbStates;
  }

 
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

  private:
    
    bool m_extrapolate;
  
    CFuint m_nbStates;

    /// Flag that defines if source sockets must be created instead of sink sockets 
    bool m_firstTimeCreation;
    
  protected:
  
    /// Dynamic data sockets
    Framework::DynamicDataSocketSet<>                  m_sockets;
    Framework::DynamicDataSocketSet<Framework::GLOBAL> m_globalSockets;
    
    /// socket data handle
    Framework::DataHandle<CFreal> m_socket;
    
    
    std::string m_socketName;

}; // end of class TurbulenceFunction

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace LESDataProcessing

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_TurbulenceFunction_hh
