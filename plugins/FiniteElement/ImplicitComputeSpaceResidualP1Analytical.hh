#ifndef COOLFluiD_Numerics_FiniteElement_ImplicitComputeSpaceResidualP1Analytical_hh
#define COOLFluiD_Numerics_FiniteElement_ImplicitComputeSpaceResidualP1Analytical_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a FEM Command to be sent to Domain to be execute
 * a ComputeSpaceResidual in an implicit manner. Works only for Linear Tetrahedra, scalar problem
 */
class ImplicitComputeSpaceResidualP1Analytical : public FiniteElementMethodCom {
public:

  /// Constructor.
  explicit ImplicitComputeSpaceResidualP1Analytical(const std::string& name);

  /// Destructor.
  ~ImplicitComputeSpaceResidualP1Analytical();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /// Execute Processing actions
  void executeOnTrs();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // methods

  /// Get the nodes coordinates for the current GeometricEntity
  void getNodes(
    Framework::GeometricEntity& cell,
    CFreal* _x, CFreal* _y, CFreal* _z);

  /// Get the states for the current GeometricEntity
  void getStates(
    Framework::GeometricEntity& cell,
    RealVector& _states, CFuint* _statesIDs, bool* _statesIPU );

  /// Calculate volume of Tetrahedra given nodal coordinates
  void calcVolume(
    const CFreal* x, const CFreal* y, const CFreal* z,
    CFreal& _V);

  /// Calculate Tetrahedra contributions
  void calcTetrahedraKeFe(
    const CFreal V, const CFreal* b, const CFreal* c, const CFreal* d,
    RealMatrix& Ke, RealVector& Fe);

  /// Calculate Tetrahedra shape functions coefficients
  void calcTetrahedraShapeFunctions(
    const CFreal& V, const CFreal* x, const CFreal* y, const CFreal* z,
    CFreal* _b, CFreal* _c, CFreal* _d);


private: // data

  /// Value of conductivity
  CFreal _conductivity;

protected: // data

  /// socket for Rhs
  Framework::DataSocketSink<
                              CFreal> socket_rhs;

}; // class ImplicitComputeSpaceResidualP1Analytical

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ImplicitComputeSpaceResidualP1Analytical_hh

