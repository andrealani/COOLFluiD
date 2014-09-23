#ifndef COOLFluiD_FluctSplit_MixLayer3DRandomPerturbation_hh
#define COOLFluiD_FluctSplit_MixLayer3DRandomPerturbation_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SubSystemStatus.hh"
#include "Framework/DataProcessingData.hh"
//#include "Environment/FileHandlerInput.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class VarSetTransformer;
  }
  
    namespace Physics {
      
      namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

class MixLayer3DRandomPerturbation : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  MixLayer3DRandomPerturbation(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MixLayer3DRandomPerturbation();

 /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();


  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  
  
  

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

protected: // functions

  /// Execute this command on the TRS
  void execute();

private: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// solution to update variable transformation
  Common::SelfRegistPtr<Framework::VarSetTransformer> _solToUpdateVar;
  
  /// update to solution variable transformation
  Common::SelfRegistPtr<Framework::VarSetTransformer> _updateToSolVar;
  
  CFreal m_Velocity; //Freestream velocity of symmetric layer
  CFreal m_di; // initial vorticity thickness
  CFreal m_order; //order relative to velocity 
  
}; 

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif 
