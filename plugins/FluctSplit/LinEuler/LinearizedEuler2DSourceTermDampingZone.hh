#ifndef COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermDampingZone_hh
#define COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermDampingZone_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"
#include "FluctSplitLinEuler.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;
  }

  namespace Physics
  {
    namespace  LinearizedEuler {   class LinEuler2DVarSet;  }
  }


    namespace FluctSplit {

      class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a source term for the LinearizedEuler2D damping zone
/// @author Erik Torres
class FluctSplitLinEuler_API LinearizedEuler2DSourceTermDampingZone : public ComputeSourceTermFSM {
public:

  /// Constructor
  /// @see LinEuler2D
  LinearizedEuler2DSourceTermDampingZone(const std::string& name);

  /// Default destructor
  ~LinearizedEuler2DSourceTermDampingZone();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Compute the source term
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
        RealVector& source,
        const FluctSplit::InwardNormalsData& normalsData);

  /// Returns the DataSockets that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


private: // data

  /// the socket stores the data of the damping coefficient
  Framework::DataSocketSink<CFreal> socket_dampingCoeff;

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// vector to store temporary result
  RealVector _temp;

  /// corresponding variable set
  Common::SafePtr<Physics::LinearizedEuler::LinEuler2DVarSet> _varSet;

}; // end of class LinEuler2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace FluctSplit

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermDampingZone_hh
