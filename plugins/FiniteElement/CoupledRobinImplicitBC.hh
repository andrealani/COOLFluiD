#ifndef COOLFluiD_Numerics_FiniteElement_CoupledRobinImplicitBC_hh
#define COOLFluiD_Numerics_FiniteElement_CoupledRobinImplicitBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElement/CoupledNeumannImplicitBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Robin Boundary condition command
 * for use when coupling subsystems
 *
 * @author Thomas Wuilbaut
 *
 */
class CoupledRobinImplicitBC : public CoupledNeumannImplicitBC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledRobinImplicitBC(const std::string& name);

  /**
   * Default destructor
   */
  ~CoupledRobinImplicitBC();

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


protected: // data

  /// a vector of string to hold the functions
  std::vector<std::string> _functionsRobin;

  /// a vector of string to hold the functions
  std::vector<std::string> _varsRobin;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunctionRobin;

  /// type of Robin BC to use
  std::string _robinType;

}; // end of class CoupledRobinImplicitBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_CoupledRobinImplicitBC_hh
