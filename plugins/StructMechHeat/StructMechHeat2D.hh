#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeat2D_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeat2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "StructMechHeatPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
  }

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a StructMechHeat2D.
 *
 * @author Thomas Wuilbaut
 *
 */
class StructMechHeat2D : public StructMechHeatPhysicalModel {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  StructMechHeat2D(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMechHeat2D();

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
   * @return the local thickness of the body
   */
  CFreal getThickness(const RealVector& coord);

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Gets the axis of symmetry
   */
  bool isAxisymmetric()
  {
    return m_isAxisymmetric;
  }

  /**
   * Gets the axis of symmetry
   */
  std::string getAxiSymmetryAxis()
  {
    return m_axisymmetryAxis;
  }

private:

  //Thickness
  CFreal m_thickness;

  //Flag if axisymmetric
  bool m_isAxisymmetric;

  //Name of the symmetry axis (X or Y)
  std::string m_axisymmetryAxis;

  //Index of the coordinate that is the radius
  CFuint m_radiusID;

}; // end of class StructMechHeat2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMechHeat_StructMechHeat2D_hh
