#ifndef COOLFluiD_Physics_Chemistry_CH4_ChemCH43D_hh
#define COOLFluiD_Physics_Chemistry_CH4_ChemCH43D_hh

//////////////////////////////////////////////////////////////////////////////

#include "ChemCH4PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the implementation of a ChemCH4 3D Physical Mocel
 *
 * @author Tiago Quintino
 *
 */
class ChemCH43D : public ChemCH4PhysicalModel {
public:

  /**
   * Constructor without arguments
   */
  ChemCH43D(const std::string& name);

  /**
   * Default destructor
   */
  ~ChemCH43D();

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * Check if this state is in a valid state
   */
  bool validate(const Framework::State& state) const;

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

}; // end of class ChemCH43D

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Chemistry_ChemCH43D_hh
