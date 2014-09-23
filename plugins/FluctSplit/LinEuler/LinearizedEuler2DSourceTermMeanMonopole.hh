#ifndef COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermMeanMonopole_hh
#define COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermMeanMonopole_hh

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
 * This class represents a source term for LinearizedEuler2D MONOPOLE with nonuniform flow
 *
 * @author Lilla Koloszar
 */
class LinearizedEuler2DSourceTermMeanMonopole : public ComputeSourceTermFSM {
public:

  /**
   * Constructor
   * @see LinEuler2D
   */
  LinearizedEuler2DSourceTermMeanMonopole(const std::string& name);

  /**
   * Default destructor
   */
  ~LinearizedEuler2DSourceTermMeanMonopole();

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
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

/**
  * Define configuration options in the CFcase config file
  */

  static void defineConfigOptions(Config::OptionList& options);


private: // data

  /// corresponding variable set
  Common::SafePtr<Physics::LinearizedEuler::LinEuler2DVarSet> _varSet;

  /// acquaintance of the PhysicalModel
  Common::SafePtr<Physics::LinearizedEuler::LinEulerTerm> _model;

  /// vector to store temporary result
  RealVector _temp;

  /// the socket stores the data of the mean flow
  Framework::DataSocketSink<RealVector> socket_meanflow;

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// width of the source defined interactively
  CFreal m_alpha;

  /// amplitude of the source defined interactively
  CFreal m_eps;

  /// frequency of the source defined interactively
  CFreal m_freq;

  /// x coordinate of the source location
  CFreal m_xs;

  /// y coordinate of the source location
  CFreal m_ys;

 }; // end of class LinEuler2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermMeanMonopole_hh
