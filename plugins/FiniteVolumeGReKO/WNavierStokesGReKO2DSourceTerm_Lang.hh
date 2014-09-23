#ifndef COOLFluiD_Numerics_FiniteVolume_WNavierStokesGReKO2DSourceTerm_Lang_hh
#define COOLFluiD_Numerics_FiniteVolume_WNavierStokesGReKO2DSourceTerm_Lang_hh

//////////////////////////////////////////////////////////////////////////////
#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Physics {

    namespace GReKO {
      class NavierStokes2DGReKOPuvt;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the 2D Gamma-Re_theta-K-Omega Source Term
 *
 * @author Khalil Bensassi 
 *
 */
class WNavierStokesGReKO2DSourceTerm_Lang : public ComputeSourceTermFVMCC {

public:

  typedef Physics::NavierStokes::MultiScalarVarSet
  <Physics::NavierStokes::Euler2DVarSet> Euler2DGReKOVarSet;
    
    /**                                                                      
    * Defines the Config Option's of this class                               
    * @param options a OptionList where to add the Option's                   
    */                                                                        



  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  WNavierStokesGReKO2DSourceTerm_Lang(const std::string& name);

  /**
   * Default destructor
   */
  ~WNavierStokesGReKO2DSourceTerm_Lang();
 


  //static void COOLFluiD::Numerics::FiniteVolume::WNavierStokesGReKO2DSourceTerm_Lang::defineConfigOptions(COOLFluiD::Config::OptionList&);
  static void defineConfigOptions(Config::OptionList& options);


  /**
   * Configure the object
   */
  virtual void configure (Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);

    _sockets.createSocketSink<RealVector>("nstates");
    _sockets.createSocketSink<CFreal>("wallDistance");
  }

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();
 
  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);
  
private: // data
  /// corresponding variable set
  Common::SafePtr<Euler2DGReKOVarSet> _varSet;

  /// corresponding diffusive variable set
  Common::SafePtr<Physics::GReKO::NavierStokes2DGReKOPuvt> _diffVarSet;

  /// vector to store temporary result
  RealVector _temp;

  /// average State
  RealVector _avState;

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

  ///K the Farfield
  CFreal _kamb;

  ///Omega at the Farfield
  CFreal _omegaamb;
  
 // Flag is true if we are using Vortivity source term
   // instead of the exact source term
  bool _SST_V;  
   // Flag is true if we are adding a sustaining terma in order
   //  to control the non-physical decay of turbulence variables 
   //  in the freestream
  bool _SST_sust;  

  RealVector   _Theta;
  RealVector   _Meanalue;

 CFreal  _Flambdak; 
 CFreal  _Flambdakprime; 
 CFreal   _Rethetat; 
 CFreal   _Rethetac; 
 CFreal   _Flength; 
 CFuint _PGrad; 

}; // end of class WNavierStokesGReKO2DSourceTerm_Lang

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NavierStokesGReKO2DSourceTerm_hh
