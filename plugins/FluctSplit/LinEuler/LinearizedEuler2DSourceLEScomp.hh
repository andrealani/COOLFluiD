#ifndef COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceLEScomp_hh
#define COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceLEScomp_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "LinEuler/LinEuler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State; 
  }

  namespace FluctSplit {

    class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term for LinearizedEuler2D LEScomp for compressible LES simulations
 * 
 * @author  Lilla Koloszar
 */
class LinearizedEuler2DSourceLEScomp : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see LinEuler2D
   */
  LinearizedEuler2DSourceLEScomp(const std::string& name);
  
  /**
   * Default destructor
   */
  ~LinearizedEuler2DSourceLEScomp();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();


  /**
   * Read the source data
   */
  virtual void readOFSource();

  /**
   * Compute the source term
   */ 
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
				RealVector& source,
				const FluctSplit::InwardNormalsData& normalsData);

/**
  * Define configuration options in the CFcase config file
  */

  static void defineConfigOptions(Config::OptionList& options);

  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets(); 
  //sourceData;

  
private: // data
  
  /// corresponding variable set
  Common::SafePtr<Physics::LinearizedEuler::LinEuler2DVarSet> _varSet;

  /// acquaintance of the PhysicalModel
  Common::SafePtr<Physics::LinearizedEuler::LinEulerTerm> _model;
  
  /// vector to store temporary result
  RealVector _temp;

  //generic input file name
  std::string m_InputFile;
  std::string m_prePath;
  std::string m_postPath;
  std::string m_postMeanPath;
  CFreal m_StartTime;
  CFreal m_TimeStep;
  CFreal m_radius;
  CFuint m_nbIter;
  bool m_ComputeSources;
  bool m_Interpolation;
  bool m_CTable;
  bool readMean_flag;


  ///mean density
  CFreal m_meanDensity;  
  CFuint nbstatesLEScomp;
  std::vector< int > CTable;
  
  ///source reading storage
  //std::vector< RealVector > Si;
  //std::vector< RealVector > SiLEE;
  
  ///past time
  CFreal m_pastTime;

 /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;
  
    /// the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;
  
    /// the socket stores the data of the velocity perturbations
  std::vector< RealVector > sourceData;
  
    /// the socket stores the data of the velocity perturbations
  Framework::DataSocketSource<RealVector> socket_sourceSave;

 }; // end of class LLinearizedEuler2DSourceLEScomp

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceLEScomp_hh
