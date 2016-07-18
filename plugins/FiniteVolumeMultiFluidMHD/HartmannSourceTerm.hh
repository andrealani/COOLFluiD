#ifndef COOLFluiD_Numerics_FiniteVolume_HartmannSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_HartmannSourceTerm_hh

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
class HartmannSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  HartmannSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~HartmannSourceTerm();

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
    
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;
  
  /// handle to outward normal
  Framework::DataHandle<CFint> _isOutward;
  
  /// socket for storing the divergence of magnetic field
  Framework::DataSocketSource<CFreal> socket_divB;
  
  /// socket for storing the Current
  Framework::DataSocketSource<CFreal> socket_Current;

  /// socket for storing the Bx potential
  Framework::DataSocketSource<CFreal> socket_BxPotential;
 
  /// socket for storing the By potential
  Framework::DataSocketSource<CFreal> socket_ByPotential;

  /// socket for storing the Bz potential
  Framework::DataSocketSource<CFreal> socket_BzPotential;

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
  
private:

  /// Electrical conductivity
  CFreal _electricalConductivity;

  /// Orszag-Tang Conductivity Flag
  bool _isResistive;

  /// Left State to compute the Orszag Conductivity
  RealVector _dataLeftState;

   /// Right State to compute the Orszag Conductivity
  RealVector _dataRightState;

  /// gradient of Bx
  RealVector _gradBx;

  /// gradient of By
  RealVector _gradBy;

 /// gradient of Bz
 RealVector _gradBz;

}; // end of class HartmannSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "HartmannSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HartmannSourceTerm_hh
