#ifndef COOLFluiD_Physics_Heat_Heat2D_hh
#define COOLFluiD_Physics_Heat_Heat2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "HeatPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a Heat2D.
 *
 * @author Tiago Quintino
 *
 */
class Heat2D : public HeatPhysicalModel {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  Heat2D(const std::string& name);

  /**
   * Default destructor
   */
  ~Heat2D();

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * @return the local thickness of the body
   */
  CFreal getThickness(const RealVector& coord);

  /**
   * Check if this state is in a valid state
   */
  bool validate(const Framework::State& state) const;

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

private:

  //Thickness
  CFreal m_thickness;

  //Flag if axisymmetric
  bool m_isAxisymmetric;

  //Name of the symmetry axis (X or Y)
  std::string m_axisymmetryAxis;

  //Index of the coordinate that is the radius
  CFuint m_radiusID;

}; // end of class Heat2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat2D_hh
