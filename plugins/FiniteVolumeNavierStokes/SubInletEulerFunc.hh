#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEulerFunc_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEulerFunc_hh

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
   * @author Radek Honzatko
   *
   */
class SubInletEulerFunc : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEulerFunc(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubInletEulerFunc();

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
  virtual void setGhostState(Framework::GeometricEntity *const face) = 0;
  
protected: // functions
  
  /// @return the area of the injector
  CFreal getArea() const
  {
    if (_inletRadii.size() == 2) {
      return MathTools::MathConsts::CFrealPi()*(_inletRadii[1]*_inletRadii[1] - _inletRadii[0]*_inletRadii[0]);
    }
    cf_assert(_inletRadii.size() == 1 && _width > 0. && _nbRings > 0);
    return 2.* MathTools::MathConsts::CFrealPi()*_inletRadii[0]*_width*(CFreal)_nbRings;
  }
  
protected: // data
  
  /// physical model data
  RealVector                               _dataInnerState;

  /// checks if an function is used in the inlet
  bool                                   _useFunction;

  /// array for temporary mass flow, v, T
  RealVector                               _inletData;

  /// storage for the temporary boundary point coordinates
  RealVector                               _bCoord;
  
  /// mass flow
  CFreal                                   _massFlow;

  /// static temperature
  CFreal                                   _temperature;
  
  /// width of the inlet injector surface in 3D
  CFreal                                   _width;
  
  /// number of rings of the inlet injector surface in 3D
  CFuint                                   _nbRings;
  
  /// a vector of radii of the torch inlet
  std::vector<CFreal>                      _inletRadii;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _functions;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction             _vFunction;
  
  /// flag telling if to inject radially (in 3D)
  bool             _radialInjection;

}; // end of class SubInletEulerFunc

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEulerFunc_hh
