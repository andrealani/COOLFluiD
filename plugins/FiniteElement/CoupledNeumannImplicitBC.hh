#ifndef COOLFluiD_Numerics_FiniteElement_CoupledNeumannImplicitBC_hh
#define COOLFluiD_Numerics_FiniteElement_CoupledNeumannImplicitBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Framework/Storage.hh"
#include "Common/CFMap.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "CoupledNeumannEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Neumann Boundary condition command
 * for use when coupling subsystems
 *
 * @author Thomas Wuilbaut
 *
 */
class CoupledNeumannImplicitBC : public FiniteElementMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledNeumannImplicitBC(const std::string& name);

  /**
   * Default destructor
   */
  ~CoupledNeumannImplicitBC();

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
  Framework::DataSocketSink< CFreal > socket_rhs;

  /// Map with pointers to the vectors of Normals
  Common::CFMap<std::string,std::vector<RealVector>*> _mapNormals;

  ///Pointer to the LocalElementData
  Common::SafePtr<LocalElementData> _localElemData;

  ///Pointer to the NeumannEntity
  Common::SelfRegistPtr<CoupledNeumannEntity> _coupledNeumannEntity;

  /// result of the integration, unperturbed and perturbed and its derivative
  RealVector _intNormal;
  RealVector _intPertbd;
  RealVector _dInt;

  /// result of the integration
  std::string _interfaceName;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  ///Do you want to alternate with anotherBC
  bool _alternateBC;

  ///If you alternate, do you want to apply this BC at start
  bool _alternateStart;

  ///If you alternate, should you apply this BC now
  bool _currentAlternateRun;

  ///Use the Robin type of BC
  bool _isRobinBC;

  /// stored configuration arguments
  /// @todo this should be avoided and removed.
  ///       It is currently only a quick fix for delayed configuration of an object (MeshFormatConverter)
  Config::ConfigArgs m_stored_args;

}; // end of class CoupledNeumannImplicitBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_CoupledNeumannImplicitBC_hh
