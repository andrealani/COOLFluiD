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

  /// a vector of radii of the torch inlet
  std::vector<CFreal>                      _inletRadii;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _functions;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction             _vFunction;
}; // end of class SubInletEulerFunc

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEulerFunc_hh
