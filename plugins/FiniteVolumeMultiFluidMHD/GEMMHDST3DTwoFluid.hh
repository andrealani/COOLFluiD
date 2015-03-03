//***************************************************************************
//
// Title: Template to create a Source term in multi-fluid MHD module
//
// Description: This class is composed of three files:
//   1.- GEMMHDST3DTwoFluid.hh: Is the header. Only declaration of the methods (functions) and members (variables used inside the class).
//   2.- GEMMHDST3DTwoFluid.cxx: Only implemented the provider of the class. It is used to be refered in the .CFcase file.
//   3.- GEMMHDST3DTwoFluid.ci: Is the main implementation in here.
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
//     8 -- massConservation_electron
//     9 -- massConservation_ion
//     10 -- massConservation_neutral
//     11 -- momConservation_electron
//     12 -- momConservation_ion
//     13 -- momConservation_neutral
//     14 -- energyConservation_electron
//     15 -- energyConservation_ion
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
//     _physicalData[11] -- y_electrons (partial density of ions) ============= firstDensity
//     _physicalData[12] -- y_ions (partial density of electron)
//     _physicalData[13] -- y_neutral (partial density of neutral)
//     _physicalData[14] -- U_electron ======================================= firstVelocity
//     _physicalData[15] -- V_electron
//     _physicalData[16] -- U_ion
//     _physicalData[17] -- V_ion
//     _physicalData[18] -- U_neutral
//     _physicalData[19] -- V_neutral
//     _physicalData[20] -- T_electron ======================================= firstTemperature
//     _physicalData[21] -- p_electron (pressure)
//     _physicalData[22] -- a_electron (speed of sound)
//     _physicalData[23] -- H_electron (total enthalpy)
//     _physicalData[24] -- T_ion
//     _physicalData[25] -- p_ion (pressure)
//     _physicalData[26] -- a_ion (speed of sound)
//     _physicalData[27] -- H_ion (total enthalpy)
//     _physicalData[28] -- T_neutral
//     _physicalData[29] -- p_neutral (pressure)
//     _physicalData[30] -- a_neutral (speed of sound)
//     _physicalData[31] -- H_neutral (total enthalpy)
//
//
// Author: Alejandro Alvarez
// Rev History: 2/2015
//
//***************************************************************************



#ifndef COOLFluiD_Numerics_FiniteVolume_GEMMHDST3DTwoFluid_hh
#define COOLFluiD_Numerics_FiniteVolume_GEMMHDST3DTwoFluid_hh

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
class GEMMHDST3DTwoFluid : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  GEMMHDST3DTwoFluid(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~GEMMHDST3DTwoFluid();

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
  * Compute the electric Current
  */
  void computeEMField();

  
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;  
  
  /// socket for storing the Ionization Rate
  //Framework::DataSocketSource<CFreal> socket_GammaIon;
  
  /// socket for storing the Ionization Rate
  //Framework::DataSocketSource<CFreal> socket_GammaRec;  
    
  /// pointer to the physical-chemical library
  //Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
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
  
  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients; 
  
  ///Vector storing the total magnetic Field
  RealVector _Btotal;
  
  ///Vector storing the total electric Field
  RealVector _Etotal;
  
  ///rate of particles of neutrals created in ionization per unit vol
  //CFreal _GammaIon_n;
  
  ///rate of particles of ions created in recombination per unit vol
  //CFreal _GammaRec_i;
 
private:

  //options here
  
}; // end of class GEMMHDST3DTwoFluid

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "GEMMHDST3DTwoFluid.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_GEMMHDST3DTwoFluid_hh
