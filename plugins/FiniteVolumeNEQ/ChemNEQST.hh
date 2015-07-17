#ifndef COOLFluiD_Numerics_FiniteVolume_ChemNEQST_hh
#define COOLFluiD_Numerics_FiniteVolume_ChemNEQST_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }
  
  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 2D for conservative
 * variables
 *
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class ChemNEQST : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  ChemNEQST(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ChemNEQST();

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
  }
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian);
  
protected:
  
  /**
   * Set the adimensional vibrational temperatures
   */
  virtual void setVibTemperature(const RealVector& pdata, 
				 const Framework::State& state,
				 RealVector& tvib);
  /**
   * Compute the source term for the axisymmetric Navier-Stokes
   */
  void computeAxiNS(Framework::GeometricEntity *const element,
		    RealVector& source);
  
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// corresponding diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffVarSet;
  
  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;
  
  /// handle to outward normal
  Framework::DataHandle<CFint> _isOutward;  
   
  /// handle to qrad
  Framework::DataHandle<CFreal> _qrad;  
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
    
  /// array to store the production/destruction term
  RealVector _omega;
  
  /// array to store the mass fractions
  RealVector _ys;
  
  /// Euler physical data
  RealVector _physicalData;

  /// vibrational temperatures
  RealVector _tvDim;
    
  /// vector to store temporary result
  RealVector _temp;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;

  /// array of temporary values
  RealMatrix _values;

  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients;
  
  /// flag telling if to use radiation coupling
  bool _hasRadiationCoupling;
  
  /// flag telling if to compute the axisymmetric source term for Navier-Stokes
  bool _includeAxiNS;
  
  /// vector of IDs for u and v components
  std::vector<CFuint> _uvID;
}; // end of class ChemNEQST

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ChemNEQST.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ChemNEQST_hh
