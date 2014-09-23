#ifndef COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermDampingMonopole_hh
#define COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermDampingMonopole_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "LinEuler/LinEuler2DVarSet.hh"
#include "FluctSplitLinEuler.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework
  {
    class State;
  }

  namespace Physics
  {
    namespace  LinearizedEuler {   class LinEuler2DVarSet;  }
  }


  namespace FluctSplit
  {

    class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a source term for LinearizedEuler2D MONOPOLE with Damping zone
/// @author  Lilla Edit Koloszar
/// @author  Erik Torres
class FluctSplitLinEuler_API LinearizedEuler2DSourceTermDampingMonopole : public ComputeSourceTermFSM {
  public:

  /// Constructor
  /// @see LinEuler2D
  LinearizedEuler2DSourceTermDampingMonopole(const std::string& name);

  /// Default destructor
  ~LinearizedEuler2DSourceTermDampingMonopole();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();
  
  /// Compute the source term
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
  			RealVector& source,
  			const FluctSplit::InwardNormalsData& normalsData);

  /// Define configuration options in the CFcase config file

  static void defineConfigOptions(Config::OptionList& options);

  /// Returns the DataSockets that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


private: // data

  /// the socket stores the data of the damping coefficient
  Framework::DataSocketSink<CFreal> socket_dampingCoeff;

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;
  /// corresponding variable set
  Common::SafePtr<Physics::LinearizedEuler::LinEuler2DVarSet> _varSet;

  /// acquaintance of the PhysicalModel
  Common::SafePtr<Physics::LinearizedEuler::LinEulerTerm> _model;

  /// vector to store temporary result
  RealVector _temp;

  /// width of the source
  CFreal m_alpha;

  /// amplitude of the source
  CFreal m_eps;

  /// frequency of the source
  CFreal m_freq;

  /// location vector of the source
  std::vector<CFreal> _sourceloc;

}; // end of class LinEuler2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace FluctSplit

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermDampingMonopole_hh
