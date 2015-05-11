#ifndef COOLFluiD_Numerics_FiniteVolume_FarFieldEuler2D_hh
#define COOLFluiD_Numerics_FiniteVolume_FarFieldEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "FiniteVolume/SuperInletInterp.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to implement farfield boundary condition
 * By deriving from @see SuperInletInterp, it allows to interpolate data from 
 * a given file 
 *
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 * @author Thomas Wuilbaut
 *
 */
class FarFieldEuler2D : public SuperInletInterp {
public:

  typedef Framework::VarSetTransformer VectorTransformer;

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FarFieldEuler2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FarFieldEuler2D();

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
  
protected: // data
  
  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// physical model data
  RealVector _pdata;

  /// array for temporary P,u,v,T
  RealVector                               _PuvT;

  /// checks if an function is used in the inlet
  bool                                   _useFunction;

  /// storage for the temporary boundary point coordinates
  RealVector                               _bCoord;

  /// static temperature
  CFreal _temperature;

  /// static pressure
  CFreal _pressure;

  /// x velocity
  CFreal _uInf;

  /// y velocity
  CFreal _vInf;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _functions;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction             _vFunction;
  
  /// RealVector holding the input variables
  Framework::State* _input;

  /// a string to hold the name of the update variables
  std::string _updateVarStr;
  
}; // end of class FarFieldEuler2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FarFieldEuler2D_hh
