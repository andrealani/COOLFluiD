#ifndef COOLFluiD_Numerics_FiniteElement_CoupledRobinBC_hh
#define COOLFluiD_Numerics_FiniteElement_CoupledRobinBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "CoupledNeumannBC.hh"

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
class CoupledRobinBC : public CoupledNeumannBC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledRobinBC(const std::string& name);

  /**
   * Default destructor
   */
  ~CoupledRobinBC();

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

  ///Type of Robin BC to use
  std::string _robinType;

}; // end of class CoupledRobinBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_CoupledRobinBC_hh
