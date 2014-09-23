#ifndef COOLFluiD_Numerics_FiniteElement_MeshMovementDirichletBC_hh
#define COOLFluiD_Numerics_FiniteElement_MeshMovementDirichletBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Framework/Storage.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Dirichlet Boundary condition command for
 * mesh movement algorithm
 *
 * @author Thomas Wuilbaut
 *
 */
class MeshMovementDirichletBC : public FiniteElementMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  MeshMovementDirichletBC(const std::string& name);

  /**
   * Default destructor
   */
  ~MeshMovementDirichletBC();

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

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on the current TRS
   */
  void executeOnTrs();

protected: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  /// socket for Rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for isUpdated
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for
  Framework::DataSocketSink<
                            std::valarray<Framework::State*> > socket_bStatesNeighbors;

  /// is this explicit or implicit
  bool _isImplicit;

}; // end of class MeshMovementDirichletBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_MeshMovementDirichletBC_hh
