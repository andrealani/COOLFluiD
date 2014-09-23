#ifndef COOLFluiD_FluctSplit_LoadSourceData3D_hh
#define COOLFluiD_FluctSplit_LoadSourceData3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SubSystemStatus.hh"
#include "Framework/DataProcessingData.hh"
#include "Environment/FileHandlerInput.hh"
#include "Framework/PhysicalModel.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
    namespace Physics {

      namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

class LoadSourceData3D : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  LoadSourceData3D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LoadSourceData3D();

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

  //generic input file name
  std::string m_InputFile;
  CFuint m_StartTime;
  CFreal m_TimeStep;
  CFuint m_nbIter;

  /// the socket stores the data of the velocity perturbations
  Framework::DataSocketSource<RealVector> socket_sourceData;

}; 

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler
    
   } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif 
