#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEulerPvtVnT_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEulerPvtVnT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command where the velocity is normal to the surface at any point
   *
   * @author Andrea Lani
   * @author Radek Honzatko
   * @author Janos Molnar
   *
   */
template <class VARSET>
class SubInletEulerPvtVnT : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEulerPvtVnT(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletEulerPvtVnT();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private: // data
  
  /// physical model var set
  Common::SafePtr<VARSET> _varSet;
  
  /// array for temporary u,v,w,T
  RealVector                               _uvwT;

  /// array for temporary Vn,T
  RealVector                               _VnT;
  
  /// array for reference values for adimensionalization of u,v,w,T
  RealVector                               _uvwTRef;

  /// storage for the temporary boundary point coordinates
  RealVector                               _bCoord;
  
  /// a vector of string to hold the functions
  std::vector<std::string>                    _functions;

  /// a vector of string to hold the functions
  std::vector<std::string>                    _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction             _vFunction;

}; // end of class SubInletEulerPvtVnT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubInletEulerPvtVnT.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEulerPvtVnT_hh
