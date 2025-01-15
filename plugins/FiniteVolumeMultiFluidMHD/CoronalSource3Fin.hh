#ifndef COOLFluiD_Numerics_FiniteVolume_CoronalSource3Fin_hh
#define COOLFluiD_Numerics_FiniteVolume_CoronalSource3Fin_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"

#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Source term for MultiFluid considering 2 fluids: plasma + neutrals
 * variables
 *
 * @author Michaela Brchnelova 
 *
 */
template <class UPDATEVAR>
class CoronalSource3Fin : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  CoronalSource3Fin(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CoronalSource3Fin();

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
  
  /**
   * Returns the DataSocket's that this command provides as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Compute the mass and reaction energy source terms
   */

  void computeMassReactionsEnergySourceTerm(float rhoe, float rhoi, float rhon, float ue, float ve, float we, float ui, float vi, float wi, float un, float vn, float wn, float Te, float Ti, float Tn);



  void computeSpitzerResistivity();



  /**
  * Compute the electric Current
  */
  void computeElectricCurrent(float rhoi, float rhoe, float ui, float vi, float wi, float ue, float ve, float we); 
  
  /**
  * Compute the Collisional momentum and energy source terms
  */
  void computeCollisionalMomentumEnergy();

  /**
 *   * Compute the Collisional momentum and energy source terms
 *     */

 /*
 * YM: Compute the charge exchange momentum and energy source terms
 */
  void computeChargeExchangeTerms();

 /*
 * YM: Compute the gravity momentum and energy source terms
 */
  void computegravityTerms(float x, float y, float z, float r, float rhoe, float rhoi, float rhon, float ue, float ve, float we, float ui, float vi, float wi, float un, float vn, float wn);

  
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;
  
  /// handle to outward normal
  //Framework::DataHandle<CFint> _isOutward;  
  
  /// socket for storing the pressure gradient in y
  Framework::DataSocketSource<CFreal> socket_gradPyIon;
  
  /// socket for storing the pressure gradient in y
  Framework::DataSocketSource<CFreal> socket_gradPyNeutral;

  /// socket for storing the pressure gradient in x
  Framework::DataSocketSource<CFreal> socket_gradPxIon;

  /// socket for storing the pressure gradient in x
  Framework::DataSocketSource<CFreal> socket_gradPxNeutral;

  /// socket for storing the pressure gradient in z
  Framework::DataSocketSource<CFreal> socket_gradPzIon;

  /// socket for storing the pressure gradient in z
  Framework::DataSocketSource<CFreal> socket_gradPzNeutral;

    
  /// pointer to the physical-chemical library
  //Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
    
  /// array to store the mass fractions
  RealVector _ys;
  
  /// Euler physical data
  RealVector _physicalData;

  /// vector to store temporary result
  RealVector _temp;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;
  
  /// array of temporary values
  RealMatrix _values;

  ///Non Induced Part of the electrocmagnetic Field
  RealVector _NonInducedEMField;
  
  /// Current density vector
  RealVector _J;
  
  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients;
  
  ///Vector storing the mass Source term
  RealVector _massSource;
  
  ///Vector storing the Collisional Momentum Source term
  RealVector _collMomentumSource;  

  RealVector _collMomentumSourceElectric;
  ///YGM: Gravity Momentum Source term
  RealVector _gravityMomentumSource;  

  ///YGM: Gravity Energy Source term
  RealVector _gravityEnergySource;

  ///YGM: Vector storing the Charge Exchage Source terms in the Momentum eqn
  RealVector _chargexMomentumSource;  
  
  ///YGM: Vector storing the Charge Exchage Source terms in the Energy eqn
  RealVector _chargexEnergySource; 
  
  ///Vector storing the Collisional Energy Source term
  RealVector _collEnergySource;  
  RealVector _collEnergySourceElectric;

  ///Vector storing the Collisional Energy Source term
  RealVector _ReactEnergySource;
  
  ///Vector storing the total magnetic Field
  RealVector _Btotal;
  
  ///Vector storing the total electric Field
  RealVector _Etotal;
  
  ///rate of particles of neutrals created in ionization per unit vol
  CFreal _GammaIon_n;   
  
  ///rate of particles of ions created in recombination per unit vol
  CFreal _GammaRec_i;     
 
  ///Spitzer resistivity
  CFreal SpitzerRes; 
private:

  /// Electrical conductivity
  CFreal _electricalResistivity;

  /// Using Spitzer resistivity
  bool _isSpitzer;

  /// flag for chemical reactions
  bool _chemicalReactions;

  /// flag for collisionall terms
  bool _collisionalTerms;

  /// flag for charge exchange terms
  bool _chargeExchangeTerms;

  /// flag for gravity terms
  bool _gravityTerms;

  /// flag for radiation
  bool _radiation;

  /// flag for heating functions
  bool _heating;

  CFreal _gravityPercentage;
  CFreal _currentPercentage;
  CFreal _mptome;
  // flag to only use the ion/neutral terms that Leake uses in Leake et al. 2014
  bool _leake14;


  // flag to impose dynamics on electrons 
  bool _ignoreEM;
  bool _ignoreHall;
  bool _rhoApprox;
  // the value of the electric charge
  CFreal _qE;

  CFreal _current;
  CFreal _currentAngle;

}; // end of class CoronalSource3Fin

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CoronalSource3Fin.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CoronalSource3Fin_hh
