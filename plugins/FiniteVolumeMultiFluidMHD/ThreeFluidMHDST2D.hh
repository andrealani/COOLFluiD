//***************************************************************************
//
// Title: Template to create a Source term in multi-fluid MHD module
//
// Description: This class is composed of three files:
//   1.- ThreeFluidMHDST2D.hh: Is the header. Only declaration of the methods (functions) and members (variables used inside the class).
//   2.- ThreeFluidMHDST2D.cxx: Only implemented the provider of the class. It is used to be refered in the .CFcase file.
//   3.- ThreeFluidMHDST2D.ci: Is the main implementation in here.
//
//  This class implements a source term. In CoolFluid, we solve equations of the form:
//
//  \frac{\partial U}{\partial t} + \nabla \cdot F^c = \nabla \cdot F^d + S
//
//  where U is the array of variables, F^c is the convective fluxes, F^d is the diffusive fluxes and S the source.
//
//  In here we implement the source in the array source[] of function computeSource().
//
//  NOTE: Don't forget to multiply the source times the volume since it is a Finite Volume code
//
//
// In the following we will assume:
//
//  - The order of the equations are:
//     0 -- Bx
//     1 -- By
//     2 -- Bz
//     3 -- Ex
//     4 -- Ey
//     5 -- Ez
//     6 -- Psi
//     7 -- Phi
//     8 -- massConservation_ion
//     9 -- massConservation_electron
//     10 -- massConservation_neutral
//     11 -- momConservation_ion
//     12 -- momConservation_electron
//     13 -- momConservation_neutral
//     14 -- energyConservation_ion
//     15 -- energyConservation_electron
//     16 -- energyConservation_neutral
//
//  e.g., the source for the Bx equation is: source[0].
//
//  - The data is stored in an array called: _physicalData[] as follows:
//     _physicalData[0] -- Bx
//     _physicalData[1] -- By
//     _physicalData[2] -- Bz
//     _physicalData[3] -- Ex
//     _physicalData[4] -- Ey
//     _physicalData[5] -- Ez
//     _physicalData[6] -- Psi
//     _physicalData[7] -- Phi
//     _physicalData[8] -- RHO (total density)
//     _physicalData[9] -- x (coordinate of the state)
//     _physicalData[10] -- y (coordinate of the state)
//     _physicalData[11] -- z (coordinate of the state)
//     _physicalData[12] -- y_ion (partial density of ions) ============= firstDensity
//     _physicalData[13] -- y_electron (partial density of electron)
//     _physicalData[14] -- y_neutral (partial density of neutral)
//     _physicalData[15] -- U_ion ======================================= firstVelocity
//     _physicalData[16] -- V_ion
//     _physicalData[17] -- U_electron
//     _physicalData[18] -- V_electron
//     _physicalData[19] -- U_neutral
//     _physicalData[20] -- V_neutral
//     _physicalData[21] -- T_ion ======================================= firstTemperature
//     _physicalData[22] -- p_ion (pressure)
//     _physicalData[21] -- a_ion (speed of sound)
//     _physicalData[22] -- H_ion (total enthalpy)
//     _physicalData[21] -- T_electron
//     _physicalData[22] -- p_electron (pressure)
//     _physicalData[21] -- a_electron (speed of sound)
//     _physicalData[22] -- H_electron (total enthalpy)
//     _physicalData[21] -- T_neutral
//     _physicalData[22] -- p_neutral (pressure)
//     _physicalData[21] -- a_neutral (speed of sound)
//     _physicalData[22] -- H_neutral (total enthalpy)
//
//
// Author: Alejandro Alvarez
// Rev History: 2/2015
//
//***************************************************************************



#ifndef COOLFluiD_Numerics_FiniteVolume_ThreeFluidMHDST2D_hh
#define COOLFluiD_Numerics_FiniteVolume_ThreeFluidMHDST2D_hh

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
 * This class represents a Source term for Three-fluid model considering fluids: ions + electrons + neutrals
 * variables
 *
 * @author Alejandro Alvarez
 *
 */
template <class UPDATEVAR>
class ThreeFluidMHDST2D : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  ThreeFluidMHDST2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ThreeFluidMHDST2D();

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
  
}; // end of class ThreeFluidMHDST2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ThreeFluidMHDST2D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ThreeFluidMHDST2D_hh
