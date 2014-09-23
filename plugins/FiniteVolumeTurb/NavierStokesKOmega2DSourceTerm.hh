#ifndef COOLFluiD_Numerics_FiniteVolume_NavierStokesKOmega2DSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_NavierStokesKOmega2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"

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
 * This class represents the 2D K-Omega Source Term
 *
 * @author Thomas Wuilbaut
 * @author Khalil Bensassi
 */
template <typename DIFFVARSET>
class NavierStokesKOmega2DSourceTerm : public ComputeSourceTermFVMCC {
public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
   NavierStokesKOmega2DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesKOmega2DSourceTerm();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
   static void defineConfigOptions(Config::OptionList& options);

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);

    _sockets.template createSocketSink<RealVector>("nstates");
    _sockets.template createSocketSink<CFreal>("wallDistance");
  }
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual  void SetDiffVarset();
  
  /// setup private data
  virtual  void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian);
  
  
protected: // data
  virtual void PreparecomputeSource(Framework::GeometricEntity *const element);
  
  virtual  void computeProductionTerm(const Framework::State& avState,
				      const CFreal& PcoFactor,const CFreal& MUT,
				      CFreal& KProdTerm,
				      CFreal& OmegaProdTerm);
  
  virtual void  computeDestructionTerm(const Framework::State& avState,
				       const CFreal& DcoFactor, CFreal& K_desterm, 
				       CFreal& Omega_desterm);

  virtual CFreal GetNSSourceTerm(); 

  static std::string getModuleName(); 
  
  // This trick was used by W. Dieudonne in Euphoria (with 20.)
  void LimitProductionTerms()
  {
    if (_limitProdTerms) {
      _prodTerm_k     = std::min(10.*fabs(_destructionTerm_k), _prodTerm_k);
      _prodTerm_Omega = std::min(10.*fabs(_destructionTerm_Omega), _prodTerm_Omega);
    }
  }
  
  /// corresponding variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _varSet;
  
  /// corresponding diffusive variable set
  Common::SafePtr<DIFFVARSET> _diffVarSet;
  
  /// vector to store temporary result
  RealVector _temp;

  /// Euler physical data
  RealVector _physicalData;

  /// handle to the reconstructed nodal states
  Framework::DataHandle< RealVector> _nstates;

  /// handle to the wall distance
  Framework::DataHandle< CFreal> _wallDistance;

  /// array of temporary values
  RealMatrix _values;

  /// array of temporary nodal states
  std::vector<RealVector*> _states;

  /// unperturbed Positive Part
  RealVector _unperturbedPositivePart;

  /// unperturbed Negative Part
  RealVector _unperturbedNegativePart;

  ///Vector for the gradients
  std::vector<RealVector*> _gradients;
  
  CFreal  _prodTerm_k;
  CFreal  _prodTerm_Omega;
  CFreal  _destructionTerm_Omega;
  CFreal  _destructionTerm_k;
  CFreal  _Radius;
  CFreal _volumes_elemID;
  CFreal  _vOverRadius;
  CFreal _avDist; 
  bool   _isAxisymmetric;
  
  /// W. Dieudonne's trick: limit production terms 
  bool _limitProdTerms;

}; // end of class NavierStokesKOmega2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesKOmega2DSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NavierStokesKOmega2DSourceTerm_hh
