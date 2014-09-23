#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEulerPvtVT_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEulerPvtVT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Common/BadValueException.hh"
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
template <class VARSET>
class SubInletEulerPvtVT : public FVMCC_BC {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  SubInletEulerPvtVT(const std::string& name);

  /// Default destructor
  virtual ~SubInletEulerPvtVT();

  /// Configures this object with user parameters
  virtual void configure ( Config::ConfigArgs& args );

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup ();

  /// Apply boundary condition on the given face
  void setGhostState (Framework::GeometricEntity *const face);

 private: // data
  
  /// physical model var set
  Common::SafePtr<VARSET> m_varSet;
  
  /// array for temporary u,v,w,T
  RealVector                               m_uvwT;
  
  /// array for reference values for adimensionalization of u,v,w,T
  RealVector                               m_uvwTRef;

  /// storage for the temporary boundary point coordinates
  RealVector                               m_bCoord;
  
  /// a vector of string to hold the functions
  std::vector<std::string>                    m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string>                    m_vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction             m_vFunction;

}; // end of class SubInletEulerPvtVT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SubInletEulerPvtVT.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEulerPvtVT_hh
