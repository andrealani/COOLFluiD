#ifndef COOLFluiD_Numerics_FiniteElement_NeumannBCImplicit_hh
#define COOLFluiD_Numerics_FiniteElement_NeumannBCImplicit_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Framework/Storage.hh"
#include "Common/CFMap.hh"
#include "FElemTypeData.hh"
#include "FiniteElementMethodData.hh"
#include "Framework/DataSocketSink.hh"
#include "NeumannEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Neumann Boundary condition command, implicit implementation
 *
 * @author Thomas Wuilbaut
 * @author Tiago Quintino
 * @author Geoffrey Deliege
 *
 */
class NeumannBCImplicit : public FiniteElementMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NeumannBCImplicit(const std::string& name);

  /**
   * Default destructor
   */
  ~NeumannBCImplicit();

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

protected:

  /**
   * Execute on the current TRS
   */
  void executeOnTrs();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  Common::SafePtr<LocalElementData> _localElemData;

  /// Neumann Entity
  Common::SelfRegistPtr<NeumannEntity> _neumannEntity;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  /// socket for Rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// result of the integration, unperturbed and perturbed and its derivative
  RealVector _intNormal;
  RealVector _intPertbd;
  RealVector _dInt;

  /// stored configuration arguments
  /// @todo this should be avoided and removed.
  ///       It is currently only a quick fix for delayed configuration of an object (MeshFormatConverter)
  Config::ConfigArgs m_stored_args;

}; // end of class NeumannBCImplicit

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NeumannBCImplicit_hh
