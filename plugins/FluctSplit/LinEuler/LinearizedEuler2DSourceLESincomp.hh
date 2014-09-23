#ifndef COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceLESincomp_hh
#define COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceLESincomp_hh

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
 * This class represents a source term for LinearizedEuler2D LESincomp
 * 
 * @author  Stefan Probsting
 */
class LinearizedEuler2DSourceLESincomp : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see LinEuler2D
   */
  LinearizedEuler2DSourceLESincomp(const std::string& name);
  
  /**
   * Default destructor
   */
  ~LinearizedEuler2DSourceLESincomp();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();

  /**
  * read LESincomp data from file
   */
  void readData();

  /**
   * Set the variable set
   * @pre the input pointer is non const to allow dynamic_cast
   */
  void setVarSet(Common::SafePtr<Framework::ConvectiveVarSet> varSet);
  
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
  
private: // data
  
  /// corresponding variable set
  Common::SafePtr<Physics::LinearizedEuler::LinEuler2DVarSet> _varSet;

  /// acquaintance of the PhysicalModel
  Common::SafePtr<Physics::LinearizedEuler::LinEulerTerm> _model;
  
  /// vector to store temporary result
  RealVector _temp;

  //generic input file name
  std::string m_InputFile;
  CFuint m_StartTime;
  CFreal m_TimeStep;
  CFuint m_nbIter;
  bool m_ComputeSources;
  bool m_Interpolation;

  ///mean density
  CFreal m_meanDensity;  

 /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;
  
    /// the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;
  
    /// the socket stores the data of the velocity perturbations
  Framework::DataSocketSink<RealVector> socket_sourceData;


 }; // end of class LLinearizedEuler2DSourceLESincomp

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceLESincomp_hh
