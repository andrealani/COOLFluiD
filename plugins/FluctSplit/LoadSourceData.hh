#ifndef COOLFluiD_FluctSplit_LoadSourceData_hh
#define COOLFluiD_FluctSplit_LoadSourceData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SubSystemStatus.hh"
#include "Framework/DataProcessingData.hh"
#include "Environment/FileHandlerInput.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Physics {

      namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

class LoadSourceData : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  LoadSourceData(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LoadSourceData();

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


  //generic input file name
  std::string m_InputFile;
  //number of iterations used for averaging
  CFuint m_AveragingPeriod;

  std::string m_InputType;
  CFreal m_StartTime;
CFreal m_avgStartTime;
  CFreal m_TimeStep;
  CFuint m_Interpolation;

  static RealVector m_indexTable;

 /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket stores the data of the velocity perturbations
  Framework::DataSocketSource<RealVector> socket_sourceData;

 /// the socket stores the states of the mean flow from tecplot file
  Framework::DataSocketSource<RealVector> socket_tecplot_states;

  /// the socket stores the coordinates of the mean flow from tecplot file
  Framework::DataSocketSource<RealVector> socket_tecplot_coords;


}; 

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif 
