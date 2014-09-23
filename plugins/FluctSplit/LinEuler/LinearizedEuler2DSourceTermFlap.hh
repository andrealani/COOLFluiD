#ifndef COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermFlap_hh
#define COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermFlap_hh

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
 * This class represents a source term for LinearizedEuler2D Flap VALIANT
 * 
 * @author  Lilla Edit Koloszar
 */
class LinearizedEuler2DSourceTermFlap : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see LinEuler2D
   */
  LinearizedEuler2DSourceTermFlap(const std::string& name);
  
  /**
   * Default destructor
   */
  ~LinearizedEuler2DSourceTermFlap();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();

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
  
  
 /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  
private: // data
  
  /// corresponding variable set
  Common::SafePtr<Physics::LinearizedEuler::LinEuler2DVarSet> _varSet;

  /// acquaintance of the PhysicalModel
  Common::SafePtr<Physics::LinearizedEuler::LinEulerTerm> _model;
  
  /// the socket stores the data of the mean flow
  Framework::DataSocketSink<RealVector> socket_meanflow;
  
  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;  
  
  /// vector to store temporary result
  RealVector _temp;

  /// width of the source
  CFreal m_alpha1;

  /// amplitude of the source
  CFreal m_eps1;

  /// frequency of the source
  CFreal m_freq1;

 /// x coordinate of the source location
  CFreal m_xs1;

  /// y coordinate of the source location
  CFreal m_ys1;

  /// width of the source
  CFreal m_alpha2;

  /// amplitude of the source
  CFreal m_eps2;

  /// frequency of the source
  CFreal m_freq2;

 /// x coordinate of the source location
  CFreal m_xs2;

  /// y coordinate of the source location
  CFreal m_ys2;

 }; // end of class LinEuler2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermFlap_hh
