#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVT_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command with the initial conditions given for tTotal, pTotal and alpha
   *
   * @author Andrea Lani
   * @author Radek Honzatko
   *
   */
template <class UPDATEVAR>
class SubInletEuler2DPuvtUVT : public FVMCC_BC {
public:
  
  typedef UPDATEVAR VARSET;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEuler2DPuvtUVT(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubInletEuler2DPuvtUVT();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);
  
 private: // data

  /// checks if an function is used in the inlet
  bool                                   _useFunction;

  /// array for temporary u,v,T
  RealVector                               _uvT;

  /// array for reference values for adimensionalization of u,v,T
  RealVector                               _uvTRef;

  /// storage for the temporary boundary point coordinates
  RealVector                               _bCoord;

  /// x velocity
  CFreal                                   _uinf;

  /// y velocity
  CFreal                                   _vinf;

  /// static temperature
  CFreal                                   _temperature;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _functions;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction             _vFunction;

}; // end of class SubInletEuler2DPuvtUVT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubInletEuler2DPuvtUVT.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVT_hh
