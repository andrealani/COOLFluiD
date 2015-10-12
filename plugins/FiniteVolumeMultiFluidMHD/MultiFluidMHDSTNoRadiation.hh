#ifndef COOLFluiD_Numerics_FiniteVolume_MultiFluidMHDSTNoRadiation_hh
#define COOLFluiD_Numerics_FiniteVolume_MultiFluidMHDSTNoRadiation_hh

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
 * @author Alejandro Alvarez
 *
 */
template <class UPDATEVAR>
class MultiFluidMHDSTNoRadiation : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  MultiFluidMHDSTNoRadiation(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MultiFluidMHDSTNoRadiation();

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
  void computeMassReactionsEnergySourceTerm();  
  
  /**
  * Compute the electric Current
  */
  void computeElectricCurrent(); 
  
  /**
  * Compute the Collisional momentum and energy source terms
  */
  void computeCollisionalMomentumEnergy();

  /**
 *   * Compute the Collisional momentum and energy source terms
 *     */
  void computeSpitzerResistivity();
  
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;
  
  /// handle to outward normal
  //Framework::DataHandle<CFint> _isOutward;  
  
  /// socket for storing the Ionization Rate
  Framework::DataSocketSource<CFreal> socket_GammaIon;
  
  /// socket for storing the Ionization Rate
  Framework::DataSocketSource<CFreal> socket_GammaRec;  
    
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
  
  ///Vector storing the Collisional Energy Source term
  RealVector _collEnergySource;  
  
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
  
}; // end of class MultiFluidMHDSTNoRadiation

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "MultiFluidMHDSTNoRadiation.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MultiFluidMHDSTNoRadiation_hh
