#ifndef COOLFluiD_Numerics_FiniteElement_NeumannBCImplicitP1Analytical_hh
#define COOLFluiD_Numerics_FiniteElement_NeumannBCImplicitP1Analytical_hh

#include "Framework/VectorialFunction.hh"
#include "FiniteElementMethodData.hh"
#include "Framework/DataSocketSink.hh"

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {


/**
 * This class represents a Neumann Boundary condition command,
 * implicit implementation. Does pointwise integration, with both
 * element-averaged and nodal contributions approach.
 */
class NeumannBCImplicitP1Analytical : public FiniteElementMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  NeumannBCImplicitP1Analytical(const std::string& name);

  /// Default destructor
  ~NeumannBCImplicitP1Analytical();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unsetup the private data and data of the aggregated classes
   * in this command after the  processing phase
   */
  void unsetup();

  /// Configures the command.
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// Execute on the current TRS
  void executeOnTrs();


protected: // data

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  /// socket for Rhs
  Framework::DataSocketSink<
                              CFreal> socket_rhs;

private: // methods

  /// Get the nodes coordinates for the current GeometricEntity
  void getNodes(
    Framework::GeometricEntity& cell, CFuint nbNodes,
    CFreal* _x, CFreal* _y, CFreal* _z);

  /// Get the states for the current GeometricEntity
  void getStates(
    Framework::GeometricEntity& cell, CFuint nbStates,
    RealVector& _states, CFuint* _statesIDs, bool* _statesIPU );

  /// Calculate area of boundary element given the nodal coordinates
  void calcArea(
    const CFreal* x, const CFreal* y, const CFreal* z,
    CFreal& A);

private: // data

  /// flag to indicate if the contribution of the element is averaged
  bool _useAverage;

}; // end of class NeumannBCImplicitP1Analytical


    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

#endif // COOLFluiD_Numerics_FiniteElement_NeumannBCImplicitP1Analytical_hh

