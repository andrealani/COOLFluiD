#ifndef COOLFluiD_Numerics_FluctSplit_NavierStokes2DKOmegaSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_NavierStokes2DKOmegaSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class State; 
  }
  
  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
    }
  }

//   namespace Numerics {
    
    namespace FluctSplit {
      
      class InwardNormalsData;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the k-Omega source term for the Navier-Stokes
 * equations
 * 
 * @author Lilla Koloszar
 */

template <typename DIFFVARSET>
class NavierStokes2DKOmegaSourceTerm : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DKOmegaSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokes2DKOmegaSourceTerm();
  
    /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();  

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
    ComputeSourceTermFSM::configure(args);
    
    _sockets.template createSocketSink<CFreal>("wallDistance");   

  }
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual  void SetDiffVarset();  
  
  
  /**
   * Compute the source term
   */ 
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
				RealVector& source,
				const FluctSplit::InwardNormalsData& normalsData); 
  
  
//   virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
//   {
//     std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = NavierStokes2DKOmegaSourceTerm::needsSockets();
// 
//     result.push_back(&socket_wallDistance);
// 
//     return result;
//   }
  
protected: // data
  virtual void PreparecomputeSource(Framework::GeometricEntity *const cell);
  
  virtual  void computeProductionTerm(const Framework::State& avState,
				      const CFreal& PcoFactor,const CFreal& MUT,
				      CFreal& KProdTerm,
				      CFreal& OmegaProdTerm);
  
  virtual void  computeDestructionTerm(const Framework::State& avState,
				       const CFreal& DcoFactor, CFreal& K_desterm, 
				       CFreal& Omega_desterm);

  virtual CFreal GetNSSourceTerm(const Framework::State& avState); 

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

  
//   /// the socket to the data handle of the state's
//   Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
 

  /// handle to the wall distance
  Framework::DataHandle< CFreal> _wallDistance;
  
//   /// socket for the wall distance values storage
//   Framework::DataSocketSink<CFreal> socket_wallDistance; 

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
  
  // cell volume
  CFreal m_cellVolume;
  
  // average radius in the cell
  CFreal m_avRadius;
  
}; // end of class NavierStokes2DKOmegaSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  
//  } // namespace Numerics
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DKOmegaSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NavierStokes2DKOmegaSourceTerm_hh
