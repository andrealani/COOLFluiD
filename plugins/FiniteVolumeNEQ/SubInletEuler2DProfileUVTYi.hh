#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DProfileUVTYi_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DProfileUVTYi_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/Node.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      template <class BASEVS> class MultiScalarVarSet;
      class Euler2DVarSet;
    }
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
 /**
  * This class implements a subsonic inlet given velocity profile, temperature and Yi
  *
  * @author Milan Zaloudek
  *
  */
   
class SubInletEuler2DProfileUVTYi : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEuler2DProfileUVTYi(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletEuler2DProfileUVTYi();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

private: // data

  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DVarSet> > _varSet;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// physical model data
  RealVector _dataInnerState;

  /// total temperature
  CFreal     _tTotal;

  /// velocity components
  CFreal _uinf;
  CFreal _vinf;

  /// static temperature
  CFreal _temperature;

  /// Storage of coordinates
  RealVector _bCoord;

  /// Vector for U,V and T computed from the function
  RealVector _bState;

  /// Vector for coordinates
  RealVector _variables;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;
  
  /// Yi mass fraction
  std::vector<CFreal> _Yi;
//  std::vector<CFreal> _Mm;

//  RealVector _Yi;
  RealVector _Mm;
  
}; // end of class SubInletEuler2DProfileUVTYi

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DProfileUVTYiTv_hh
