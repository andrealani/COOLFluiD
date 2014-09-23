#ifndef COOLFluiD_Numerics_FluctSplit_DNS2DSourceTermDampingZone_hh
#define COOLFluiD_Numerics_FluctSplit_DNS2DSourceTermDampingZone_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;
  }

  namespace Physics
  {
    namespace  NavierStokes {   class Euler2DVarSet;  }
  }


    namespace FluctSplit {

      class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a source term for the LinearizedEuler2D damping zone
/// @author Erik Torres
class DNS2DSourceTermDampingZone : public ComputeSourceTermFSM {
public:

  /// Constructor
  /// @see LinEuler2D
  DNS2DSourceTermDampingZone(const std::string& name);

  /// Default destructor
  ~DNS2DSourceTermDampingZone();

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


  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;
  /// the socket stores the data of the damping coefficient
  Framework::DataSocketSink<CFreal> socket_dampingCoeff;

  /// vector to store temporary result
  RealVector _temp;

  /// corresponding variable set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;

  RealVector dampingCoeff;
  CFreal m_nuMax;
  CFreal m_r0;
  CFreal m_beta;
  
  ///Name of the boundary associated to the damping zone
  std::vector<std::string> m_BoundaryTRS;

}; // end of class DNS2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace FluctSplit

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_DNS2DSourceTermDampingZone_hh
