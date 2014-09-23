#ifndef COOLFluiD_FluctSplit_Vorticity_NavierStokes3Dcons_hh
#define COOLFluiD_FluctSplit_Vorticity_NavierStokes3Dcons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FileHandlerOutput.hh"
#include "Common/Trio.hh"
#include "MathTools/FunctionParser.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "NavierStokes/NavierStokes3DCons.hh"
#include "NavierStokes/Euler3DCons.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "FluctSplit/InwardNormalsData.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class samples the solution at a point for every iteration
/// @author Lilla Kapa
class Vorticity_NavierStokes3Dcons : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  Vorticity_NavierStokes3Dcons(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~Vorticity_NavierStokes3Dcons();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
 virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


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


  /**
   * Compute the values at the wall and write them to file
   */
  void computeVorticity(bool issave );

    /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileVorticity();

private: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;

  /// handle the file where sampling is dumped
  /// handle is opened in setup() closed in unsetup()
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fhandle;

  /// Index of the cell that contains the sampling point
  std::vector<CFuint> m_SamplingCellIndex;

  /// Pointer to TRS
  std::vector< Common::SafePtr< Framework::TopologicalRegionSet> > m_PointerToSampledTRS;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// socket for inward normals
  Framework::DataSocketSink< InwardNormalsData*> socket_normals;
  
 /// Output File for Wall Values
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fileVorticity;

  /// Storage for choosing when to save the wall values file
  CFuint m_saveRateVorticity;
  CFuint m_compRateVorticity;
 /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileVorticity;

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler3DCons> m_varSet;
  
  // ///flag for appending iteration
  bool m_appendIter;

  // ///flag for appending time
  bool m_appendTime;

//   /// socket for stencil
  Framework::DataSocketSink<CFreal> socket_vorticity;

}; // end of class SamplingPoint

//////////////////////////////////////////////////////////////////////////////

    }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_SamplingPoint_hh
